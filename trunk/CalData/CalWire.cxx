////////////////////////////////////////////////////////////////////////
// $Id: CalWire.cxx,v 1.36 2010/09/15  bpage Exp $
//
// CalWire class
//
// brebel@fnal.gov, echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// LArSoft includes
#include "CalData/CalWire.h"
#include "Geometry/geo.h"
#include "Filters/ChannelFilter.h"
#include "RawData/raw.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArFFT.h"


// ROOT includes
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TComplex.h>

namespace caldata{

  //-------------------------------------------------
  CalWire::CalWire(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Wire> >();
    
  }
  
  //-------------------------------------------------
  CalWire::~CalWire()
  {
  }

  //////////////////////////////////////////////////////
  void CalWire::reconfigure(fhicl::ParameterSet const& p)
  {
    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(p.get<std::string>("ResponseFile"), fResponseFile);
    fExpEndBins       = p.get< int >        ("ExponentialEndBins");
    fPostsample       = p.get< int >        ("PostsampleBins");
  }

  //-------------------------------------------------
  void CalWire::beginJob()
  {  
    
    LOG_DEBUG("CalWire") << "CalWire_plugin: Opening  Electronics Response File: " 
			 << fResponseFile.c_str();
    
    TFile f(fResponseFile.c_str());
    if( f.IsZombie() )
      mf::LogWarning("CalWire") << "Cannot open response file " 
				<< fResponseFile.c_str();
    
    TH2D *respRe       = dynamic_cast<TH2D*>(f.Get("real/RespRe")   );
    TH2D *respIm       = dynamic_cast<TH2D*>(f.Get("real/RespIm")   );
    TH1D *decayHist    = dynamic_cast<TH1D*>(f.Get("real/decayHist"));
    unsigned int wires = decayHist->GetNbinsX();
    unsigned int bins  = respRe->GetYaxis()->GetNbins();
    unsigned int bin   = 0;
    unsigned int wire  = 0;
    fDecayConstsR.resize(wires);
    fKernMapR.resize(wires);
    fKernelR.resize(respRe->GetXaxis()->GetNbins());
    const TArrayD *edges = respRe->GetXaxis()->GetXbins();
    for(int i = 0; i < respRe->GetXaxis()->GetNbins(); ++i)  {
      fKernelR[i].resize(bins);
      for(bin = 0; bin < bins; ++bin) {  
        
	const TComplex a(respRe->GetBinContent(i+1,bin+1), 
			 respIm->GetBinContent(i+1,bin+1));
	fKernelR[i][bin]=a;
      }
      for(; wire < (*edges)[i+1]; ++wire) {
	fKernMapR[wire]=i;
        fDecayConstsR[wire]=decayHist->GetBinContent(wire+1); 
      }
    }
    respRe       = dynamic_cast<TH2D*>(f.Get("sim/RespRe")   );
    respIm       = dynamic_cast<TH2D*>(f.Get("sim/RespIm")   );
    decayHist    = dynamic_cast<TH1D*>(f.Get("sim/decayHist"));
    wires = decayHist->GetNbinsX();
    bins  = respRe->GetYaxis()->GetNbins();
    fDecayConstsS.resize(wires);
    fKernMapS.resize(wires);
    fKernelS.resize(respRe->GetXaxis()->GetNbins()); 
    const TArrayD *edges1 = respRe->GetXaxis()->GetXbins();
    wire =0;
    for(int i = 0; i < respRe->GetXaxis()->GetNbins(); ++i)  {
      fKernelS[i].resize(bins);
      for(bin = 0; bin < bins; ++bin) {  
	const TComplex b(respRe->GetBinContent(i+1,bin+1), 
			 respIm->GetBinContent(i+1,bin+1));
	fKernelS[i][bin]=b;
      }
      for(; wire < (*edges1)[i+1]; ++wire) {
	fKernMapS[wire]=i;
        fDecayConstsS[wire]=decayHist->GetBinContent(wire+1); 
      }
    }
   
    f.Close();
  }

  //////////////////////////////////////////////////////
  void CalWire::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWire::produce(art::Event& evt)
  {
    
      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<double> decayConsts;  
    std::vector<int> kernMap;
    std::vector<std::vector<TComplex> > kernel; 
    //Put correct response functions and decay constants in place
    if(evt.isRealData()) {
      decayConsts=fDecayConstsR;
      kernMap=fKernMapR;
      kernel=fKernelR;
    }
    else {
      decayConsts=fDecayConstsS;
      kernMap=fKernMapS;
      kernel=fKernelS;
    }

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;

    // make a collection of Wires
    std::auto_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel(fDigitModuleLabel, digitVecHandle);

    if (!digitVecHandle->size())  return;
    mf::LogInfo("CalWire") << "CalWire:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
    
    int transformSize = fFFT->FFTSize();
    unsigned int channel(0); // channel number
    unsigned int bin(0);     // time bin loop variable
    
    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();  

    double decayConst = 0.;  // exponential decay constant of electronics shaping
    double fitAmplitude    = 0.;  //This is the seed value for the amplitude in the exponential tail fit 
    std::vector<double> holder;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    
    // loop over all wires    
    for(unsigned int rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();

      // skip bad channels
      if(!chanFilt->BadChannel(channel)) {
	holder.resize(transformSize);
	
	// uncompress the data
	raw::Uncompress(digitVec->fADC, rawadc, digitVec->Compression());
	
	for(bin = 0; bin < dataSize; ++bin) 
	  holder[bin]=(rawadc[bin]-digitVec->GetPedestal());
	// fExpEndBins only nonzero for detectors needing exponential tail fitting
	if(fExpEndBins && fabs(decayConsts[channel]) > 0.0){
	  
	  TH1D expTailData("expTailData","Tail data for fit",
			   fExpEndBins,dataSize-fExpEndBins,dataSize);
	  TF1  expFit("expFit","[0]*exp([1]*x)");
	  
	  for(bin = 0; bin < (unsigned int)fExpEndBins; ++bin) 
	    expTailData.Fill(dataSize-fExpEndBins+bin,holder[dataSize-fExpEndBins+bin]);
	  decayConst = decayConsts[channel];
	  fitAmplitude = holder[dataSize-fExpEndBins]/exp(decayConst*(dataSize-fExpEndBins));
	  expFit.FixParameter(1,decayConst);
	  expFit.SetParameter(0,fitAmplitude);
	  expTailData.Fit(&expFit,"QWN","",dataSize-fExpEndBins,dataSize);
	  expFit.SetRange(dataSize,transformSize);
	  for(bin = 0; bin < dataSize; ++bin)
	    holder[dataSize+bin]= expFit.Eval(bin+dataSize);
	}
	// This is actually deconvolution, by way of convolution with the inverted 
        // kernel.  This code assumes the response function has already been
        // been transformed and inverted.  This way a complex multiplication, rather
        // than a complex division is performed saving 2 multiplications and 
        // 2 divsions

	if( geom->DetId() == geo::kArgoNeuT )
	  {
	    fFFT->Convolute(holder,kernel[kernMap[channel]]);
	  }
	else if( geom->DetId() == geo::kMicroBooNE )
	  {
	    // Figure out which kernel to use (0=induction, 1=collection).
	    unsigned int tpc, plane, wire;
	    geom->ChannelToWire(channel, tpc, plane, wire);
	    geo::SigType_t sigtype = geom->Plane(plane, tpc).SignalType();
	    size_t k;
	    if(sigtype == geo::kInduction)
	      k = 0;
	    else if(sigtype == geo::kCollection)
	      k = 1;
	    else
	      throw cet::exception("CalWire") << "Bad signal type = " << sigtype << "\n";
	    assert(k < kernel.size());
	    
	    fFFT->Convolute(holder,kernel[k]);
	  }
	else
	      throw cet::exception("CalWire") << "Deconvolution not handled yet for this detector. " << "\n";
      } 
      
      holder.resize(dataSize,1e-5);
      //This restores the DC component to signal removed by the deconvolution.
      if(fPostsample) {
        double average=0.0;
	for(bin=0; bin < (unsigned int)fPostsample; ++bin) 
	  average+=holder[holder.size()-1-bin]/(double)fPostsample;
        for(bin = 0; bin < holder.size(); ++bin) holder[bin]-=average;
      }  
      wirecol->push_back(recob::Wire(holder,digitVec));
    }
    
    if(wirecol->size() == 0)
      mf::LogWarning("CalWire") << "No wires made for this event.";
    
    evt.put(wirecol);
    
    delete chanFilt;
    return;
  }
  
} // end namespace caldata

