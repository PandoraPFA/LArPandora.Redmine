////////////////////////////////////////////////////////////////////////
// $Id: SimWireT962.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireT962 class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// LArSoft includes
#include "DetSim/SimWireT962.h"
#include "Geometry/geo.h"
#include "Simulation/sim.h"
#include "RawData/raw.h"
#include "Utilities/DetectorProperties.h"

// ROOT includes
#include <TMath.h>
#include <TH1D.h>
#include <TFile.h>

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace detsim{

  //-------------------------------------------------
  SimWireT962::SimWireT962(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();
 
    fCompression = raw::kNone;
    std::string compression(pset.get< std::string >("CompressionType"));
    if(compression.compare("Huffman") == 0) fCompression = raw::kHuffman;    

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine(seed);
  }

  //-------------------------------------------------
  SimWireT962::~SimWireT962()
  {
    fColFieldResponse.clear();
    fIndFieldResponse.clear();
    fColShape.clear();
    fIndShape.clear();
    fChargeWork.clear();
    fElectResponse.clear();
 
    for(unsigned int i = 0; i < fNoise.size(); ++i) fNoise[i].clear();
    fNoise.clear();

  }

  //-------------------------------------------------
  void SimWireT962::reconfigure(fhicl::ParameterSet const& p) 
  {
    fResponseSet      = false;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(p.get<std::string>("ResponseFile"), fResponseFile);

    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");
    fNoiseFact        = p.get< double              >("NoiseFact");
    fNoiseWidth       = p.get< double              >("NoiseWidth");
    fLowCutoff        = p.get< double              >("LowCutoff");
    fNFieldBins       = p.get< int                 >("FieldBins");
    fCol3DCorrection  = p.get< double              >("Col3DCorrection");
    fInd3DCorrection  = p.get< double              >("Ind3DCorrection");
    fColFieldRespAmp  = p.get< double              >("ColFieldRespAmp");
    fIndFieldRespAmp  = p.get< double              >("IndFieldRespAmp");
    fShapeTimeConst   = p.get< std::vector<double> >("ShapeTimeConst");


    art::ServiceHandle<util::DetectorProperties> detprop;
    fSampleRate       = detprop->SamplingRate();
    fTriggerOffset    = detprop->TriggerOffset();

    return;
  }

  //-------------------------------------------------
  void SimWireT962::beginJob() 
  { 
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist      = tfs->make<TH1D>("Noise", ";Noise (ADC);", 1000, -10., 10.);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();
    fChargeWork.resize(fNTicks, 0.);

    fNoise.resize(100);
    for(int p = 0; p < 100; ++p){

      GenNoise(fNoise[p]);
      for(int i = 0; i < fNTicks; ++i){
	fNoiseDist->Fill(fNoise[p][i]);
      }
    }// end loop over wires

    // set field response and electronics response, then convolute them
    SetFieldResponse();
    SetElectResponse();
    ConvoluteResponseFunctions();


    return;

  }

  //-------------------------------------------------
  void SimWireT962::endJob() 
  {
  }

  //-------------------------------------------------
  void SimWireT962::produce(art::Event& evt)
  {

    //std::cout << "in SimWireT962::produce " << std::endl;

    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    unsigned int tpc        = 0;
    unsigned int plane      = 0;
    unsigned int wire       = 0;
    unsigned int signalSize = fNTicks/2;
    // vectors for working
    std::vector<short>    adcvec(signalSize, 0);	

    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(geo->Nchannels());
    for(size_t c = 0; c < chanHandle.size(); ++c){
      channels[chanHandle[c]->Channel()] = chanHandle[c];
    }
    
    // make an auto_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::auto_ptr< std::vector<raw::RawDigit>   >  digcol(new std::vector<raw::RawDigit>);
	  
    unsigned int chan = 0; 
    fChargeWork.clear();
    fChargeWork.resize(fNTicks, 0.);
	  
    art::ServiceHandle<util::LArFFT> fFFT;

    // Add all channels  
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    std::map<int,double>::iterator mapIter;      
    for(chan = 0; chan < geo->Nchannels(); chan++) {    
      
      fChargeWork.clear();    
      fChargeWork.resize(fNTicks, 0.);    

      if( channels[chan] ){      
	geo->ChannelToWire(chan, tpc, plane, wire);      

	// get the sim::SimChannel for this channel
	const sim::SimChannel* sc = channels[chan];

	// loop over the tdcs and grab the number of electrons for each
	for(size_t t = 0; t < fChargeWork.size(); ++t) 
	  fChargeWork[t] = sc->Charge(t);      

        // Convolve charge with appropriate response function 
	if(geo->Plane(plane,tpc).SignalType() == geo::kInduction)        
	  fFFT->Convolute(fChargeWork,fIndShape);      
	else fFFT->Convolute(fChargeWork,fColShape);      
      }

      // noise was already generated for each wire in the event
      // raw digit vec is already in channel order
      // pick a new "noise channel" for every channel  - this makes sure    
      // the noise has the right coherent characteristics to be on one channel   
      int noisechan = TMath::Nint(flat.fire()*(1.*(fNoise.size()-1)+0.1));
      for(unsigned int i = 0; i < signalSize; ++i){
	adcvec[i] = (short)TMath::Nint(fNoise[noisechan][i] + fChargeWork[i]);
      }

      // compress the adc vector using the desired compression scheme,
      // if raw::kNone is selected nothing happens to adcvec
      raw::Compress(adcvec, fCompression);
      
      raw::RawDigit rd(chan, signalSize, adcvec, fCompression);

      // add this digit to the collection
      digcol->push_back(rd);
    }// end loop over channels      

    evt.put(digcol);

    return;
  }

  //-------------------------------------------------
  void SimWireT962::ConvoluteResponseFunctions()
  {
    std::vector<double> col(fNTicks, 0.);
    std::vector<double> ind(fNTicks, 0.);

    unsigned int mxbin = TMath::Min(fNTicks, (int)fNElectResp + fNFieldBins);

    double sumCol = 0.;
    double sumInd = 0.;
  
    for(unsigned int i = 1; i < mxbin; ++i){
      sumCol = 0.;
      sumInd = 0.;
      for(unsigned int j = 0; j < (unsigned int)fNFieldBins; ++j){
	unsigned int k = i - j;
	if(k == 0) break;
	sumCol += fElectResponse[k]*fColFieldResponse[j];
	sumInd += fElectResponse[k]*fIndFieldResponse[j];
      }
      col[i] = sumCol;
      ind[i] = sumInd;
    
    }//end loop over bins;

    // pad out the rest of the vector with 0.  
    ind.resize(fNTicks, 0.);  
    col.resize(fNTicks, 0.);  

    // write the shapes out to a file
    art::ServiceHandle<art::TFileService> tfs;
    fColTimeShape = tfs->make<TH1D>("ConvolutedCollection",";ticks; Electronics#timesCollection",fNTicks,0,fNTicks);
    fIndTimeShape = tfs->make<TH1D>("ConvolutedInduction",";ticks; Electronics#timesInduction",fNTicks,0,fNTicks);    

    fIndShape.resize(fNTicks/2+1);  
    fColShape.resize(fNTicks/2+1);  

    // do the FFT of the shapes  
    std::vector<double> delta(fNTicks);  
    delta[0] = 1.0;  
    delta[fNTicks-1]=1.0;  

    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->AlignedSum(ind,delta,false);  
    fFFT->AlignedSum(col,delta,false);  
    fFFT->DoFFT(ind, fIndShape);  
    fFFT->DoFFT(col, fColShape);

    // check that you did the right thing  
    for(unsigned int i = 0; i < ind.size(); ++i){    
      fColTimeShape->Fill(i, col[i]);    
      fIndTimeShape->Fill(i, ind[i]);  
    }  

    fColTimeShape->Write();
    fIndTimeShape->Write();


    return;
  }

  //-------------------------------------------------
  void SimWireT962::GenNoise(std::vector<double>& noise)
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    noise.clear();
    noise.resize(fNTicks, 0.);
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);///<noise in frequency space

    double pval = 0.; 
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);
    for(int i=0; i< fNTicks/2+1; ++i){
      // exponential noise spectrum 
      pval = fNoiseFact*exp(-(double)i*binWidth/fNoiseWidth);
      // low frequency cutoff     
      lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
      // randomize 10%
      flat.fireArray(2,rnd,0,1);
      pval *= lofilter*(0.9+0.2*rnd[0]);
      // random pahse angle
      phase = rnd[1]*2.*TMath::Pi();

      TComplex tc(pval*cos(phase),pval*sin(phase));
      noiseFrequency[i] += tc;
    }
  
    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);

    noiseFrequency.clear();

    // multiply each noise value by fNTicks as the InvFFT 
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    for(unsigned int i = 0; i < noise.size(); ++i) noise[i] *= 1.*fNTicks;

    return;
  }

  //-------------------------------------------------
  void SimWireT962::SetFieldResponse()
  {

    art::ServiceHandle<geo::Geometry> geo;
 
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double xyzl[3] = {0.};
    // should always have at least 2 planes
    geo->Plane(0).LocalToWorld(xyzl, xyz1);
    geo->Plane(1).LocalToWorld(xyzl, xyz2);

    // this assumes all planes are equidistant from each other,
    // probably not a bad assumption
    double pitch = xyz2[0] - xyz1[0]; ///in cm

    fColFieldResponse.resize(fNFieldBins, 0.);
    fIndFieldResponse.resize(fNFieldBins, 0.);

    // set the response for the collection plane first
    // the first entry is 0

    // write out the response functions to the file
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    fIndFieldResp = tfs->make<TH1D>("InductionFieldResponse",";t (ns);Induction Response",fNTicks,0,fNTicks);
    fColFieldResp = tfs->make<TH1D>("CollectionFieldResponse",";t (ns);Collection Response",fNTicks,0,fNTicks);
    art::ServiceHandle<util::LArProperties> larp;
    double driftvelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature())/1000.;  
    int nbinc = TMath::Nint(fCol3DCorrection*(fabs(pitch))/(driftvelocity*fSampleRate)); ///number of bins //KP
  
    double integral = 0.;
    for(int i = 1; i < nbinc; ++i){
      fColFieldResponse[i] = fColFieldResponse[i-1] + 1.0;
      integral += fColFieldResponse[i];
    }

    for(int i = 0; i < nbinc; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
      fColFieldResp->Fill(i, fColFieldResponse[i]);
    }

    ///now the induction plane
    
    int nbini = TMath::Nint(fInd3DCorrection*(fabs(pitch))/(driftvelocity*fSampleRate));//KP
    for(int i = 0; i < nbini; ++i){
      fIndFieldResponse[i] = fIndFieldRespAmp/(1.*nbini);
      fIndFieldResponse[nbini+i] = -fIndFieldRespAmp/(1.*nbini);
    
      fIndFieldResp->Fill(i, fIndFieldResponse[i]);
      fIndFieldResp->Fill(nbini+i, fIndFieldResponse[nbini+i]);
    
    }

    fColFieldResp->Write();
    fIndFieldResp->Write();

    return;
  }

  //-------------------------------------------------
  void SimWireT962::SetElectResponse()
  {

    art::ServiceHandle<geo::Geometry> geo;
  
    fElectResponse.resize(fNTicks, 0.);
    std::vector<double> time(fNTicks,0.);

    TFile inFile(fResponseFile.c_str(),"READ");        
    if(inFile.IsZombie()){      
      mf::LogError("SimWireT962") << "Cannot open response file" << fResponseFile.c_str()
				  << " bail ungracefully";
      assert(0);    
    }    
    
    TH1D * shape = (TH1D*)inFile.Get("shape");    
    int ctr = 0;    
    unsigned int lastpos = 0;    
    double integral = 0.;    
    double holder = 0;   
    
    while(ctr < fNTicks){      
      holder = shape->GetBinContent(ctr+1);      
      fElectResponse[ctr] = holder;      
      time[ctr] = ctr*fSampleRate;      
      if(holder> 0.){	
	integral += holder;	
	lastpos = ctr;      
      }      
      
      ++ctr;    
    }//end loop over input         
    
    for(unsigned int i = 0; i < fElectResponse.size(); ++i)      
      fElectResponse[i] /= integral;
    
    fNElectResp = fElectResponse.size();

    // write the response out to a file
    art::ServiceHandle<art::TFileService> tfs;
    fElectResp = tfs->make<TH1D>("ElectronicsResponse",";t (ns);Electronics Response",fNElectResp,0,fNElectResp);
    for(unsigned int i = 0; i < fNElectResp; ++i){
      fElectResp->Fill(i, fElectResponse[i]);
    }

    fElectResp->Write();

    return;
  }

}
