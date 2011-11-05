////////////////////////////////////////////////////////////////////////
// $Id: CalWireAna.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// CalWireAna class designed to make histos.
//
//
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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "CalData/CalWireAna.h"
#include "Geometry/geo.h"
#include "Simulation/sim.h"
#include "RawData/raw.h"
#include "RecoBase/recobase.h"
#include "Filters/ChannelFilter.h"
#include "Utilities/LArFFT.h"

// ROOT includes
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

namespace caldata{

  //-------------------------------------------------
  CalWireAna::CalWireAna(fhicl::ParameterSet const& pset) : 
    fCalWireModuleLabel(pset.get< std::string >("CalWireModuleLabel")),
    fDetSimModuleLabel(pset.get< std::string >("DetSimModuleLabel"))
  {
			 
  }

  //-------------------------------------------------
  CalWireAna::~CalWireAna()
  {
  }

  //-------------------------------------------------
  void CalWireAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fDiffsR          = tfs->make<TH1F>("One timestamp diffs in RawDigits", ";#Delta ADC;",        40,   -19.5,  20.5);
    fDiffsW          = tfs->make<TH1F>("One timestamp diffs in Wires", ";#Delta ADC;",        20,   -9.5,  10.5);
    fDiffsRW         = tfs->make<TH1F>("Same timestamp diffs in RD-Wires", ";#Delta ADC;",        20,   -9.5,  10.5);
    fDiffsRWvsR      = tfs->make<TH2F>("Same timestamp diffs in RD-Wires vs R", ";#Delta ADC;",   481,-0.5,480.5,    20,   -9.5,  10.5);
    fDiffsRWgph      = tfs->make<TH1F>("Same timestamp diffs in RD-Wires gph", ";#Delta ADC;",        20,   -9.5,  10.5);
    fDiffsRWvsRgph   = tfs->make<TH2F>("Same timestamp diffs in RD-Wires vs R gph", ";#Delta ADC;",   481,-0.5,480.5,    20,   -9.5,  10.5);
    fRawSig          = tfs->make<TH2F>("One event, one channel Raw", "timestamp",        481,-0.5,480.5, 21,-0.5,  20.5);
    fWireSig         = tfs->make<TH2F>("One event, one channel Wire", "timestamp",     481,-0.5,480.5, 21, -0.5,  20.5);

    fRD_WireMeanDiff2D = tfs->make<TH2F>("Mean (Raw-CALD)-over-Raw in Window gph","Wire number",481,-0.05,480.5, 40, -1., 1.);
    fRD_WireRMSDiff2D  = tfs->make<TH2F>("RMS (Raw-CALD)-over-Raw in Window gph","Wire number",481,-0.05,480.5, 10, 0., 2.);

    fWindow       = tfs->make<TH2F>("tmax-tmin vs indMax", "ticks",  200, 0, 2000,      20,   -2.5,  60.5);
    fMin          = tfs->make<TH1F>("Value of min", "ticks",        21,   -20.5,  0.5);
    fMax          = tfs->make<TH1F>("Value of max", "ticks",        21,   0.5,  20.5);
    fRawIndPeak      = tfs->make<TH1F>("indPeakRaw", ";Induction Peaks Raw;",40,5,45);
    fRawColPeak      = tfs->make<TH1F>("colPeakRaw", ";Collection Peaks Raw;",40,5,45);
    fCalIndPeak      = tfs->make<TH1F>("indPeakCal", ";Induction Peaks Calibrated;",40,5,45);
    fCalColPeak      = tfs->make<TH1F>("colPeakCal", ";Collection Peaks Calibrated;",40,5,45);

    fIR = tfs->make<TH1F>("Raw Ind signal","time ticks",4096,0.0,4096.);
    fCR = tfs->make<TH1F>("Raw Coll signal","time ticks",4096,0.0,4096.);
    fIW = tfs->make<TH1F>("Wire Ind signal","time ticks",4096,0.0,4096.);
    fCW = tfs->make<TH1F>("Wire Coll signal","time ticks",4096,0.0,4096.);
    fNoiseHist = tfs->make<TH1F>("Noise Histogram","FFT Bins",2049,0,2049);
    fNoiseRMS = tfs->make<TH1F>("Noise RMS","RMS",25,0,2.0);

    return;

  }

  //-------------------------------------------------
  void CalWireAna::endJob() 
  {
  }

  //-------------------------------------------------
  void CalWireAna::analyze(const art::Event& evt)
  {

    // loop over the raw digits and get the adc vector for each, then compress it and uncompress it

    filter::ChannelFilter * fil = new  filter::ChannelFilter();
    art::Handle< std::vector<raw::RawDigit> > rdHandle;
    evt.getByLabel(fDetSimModuleLabel,rdHandle);
    art::Handle< std::vector<recob::Wire> > wHandle;
    evt.getByLabel(fCalWireModuleLabel,wHandle);

    art::PtrVector<recob::Wire> wvec;
    for(unsigned int i = 0; i < wHandle->size(); ++i){
      art::Ptr<recob::Wire> w(wHandle,i);
      wvec.push_back(w);
    }
    art::PtrVector<raw::RawDigit> rdvec;
    for(unsigned int i = 0; i < rdHandle->size(); ++i){
      art::Ptr<raw::RawDigit> r(rdHandle,i);
      rdvec.push_back(r);
    }

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArFFT> fft;
    unsigned int tpc           = 0;
    unsigned int plane         = 0;
    unsigned int wire          = 0;
    double pedestal = rdvec[0]->GetPedestal();
    double threshold = 9.0;
    double signalSize = rdvec[0]->Samples();
    int indChan0=64;
    int indChan1=110;
    int colChan0=312;
    int colChan1=354;
    std::vector<double> ir(fft->FFTSize()),iw(fft->FFTSize()),cr(fft->FFTSize()),cw(fft->FFTSize());
    ir[signalSize]=iw[signalSize]=cr[signalSize]=cw[signalSize]=1.0;
    /// loop over all the raw digits in a window around peak 
    for(unsigned int rd = 0; rd < rdvec.size(); ++rd){
      // Find corresponding wire.
      std::vector<double> signal(fft->FFTSize());
      for(unsigned int wd = 0; wd < wvec.size(); ++wd){
	if (wvec[wd]->RawDigit() == rdvec[rd]){
	  signal.insert(signal.begin(),wvec[wd]->fSignal.begin(),wvec[wd]->fSignal.end()); 
	  signal.resize(fft->FFTSize());
	  break;
	}
	if (wd == (wvec.size()-1) ){
	  LOG_DEBUG("CalWireAna") << "caldata::CalWireAna:Big problem! No matching Wire for RawDigit. Bailing." << rd;
	  return;
	}
      }

      geom->ChannelToWire(rdvec[rd]->Channel(), tpc, plane, wire);
      
      std::vector<double> adc(fft->FFTSize());
      
      for(unsigned int t = 1; t < rdvec[rd]->Samples(); ++t){
	fDiffsR->Fill(rdvec[rd]->ADC(t) - rdvec[rd]->ADC(t-1)); 
	adc[t-1]=rdvec[rd]->ADC(t-1);
	fRawSig->Fill(rd,rdvec[rd]->ADC(t));
      }
      //get the last one for the adc vector
      adc[rdvec[rd]->Samples()-1]=rdvec[rd]->ADC(rdvec[rd]->Samples()-1);
      if(!fil->BadChannel(rdvec[rd]->Channel()) && (*max_element(adc.begin(),adc.end()) < pedestal+threshold && *min_element(adc.begin(),adc.end()) >pedestal -threshold)) {
	double sum=0;
	for(int i = 0; i < signalSize; i++) sum+=pow(adc[i]-pedestal,2.0);
	fNoiseRMS->Fill(TMath::Sqrt(sum/(double)signalSize));
	std::vector<double> temp(fft->FFTSize());
	std::vector<TComplex> fTemp(fft->FFTSize()/2+1);
	for(int i = 0; i < signalSize; i++) temp[i]=(adc[i]-pedestal)*sin(TMath::Pi()*(double)i/signalSize);
	fft->DoFFT(temp,fTemp);
	for(int i = 0; i < fft->FFTSize()/2+1; i++) fNoiseHist->Fill(i,fTemp[i].Rho());
      }
      if(geom->Plane(plane,tpc).SignalType() == geo::kInduction && 
	 rdvec[rd]->Channel() > indChan0 && 
	 rdvec[rd]->Channel() < indChan1){	    
	    fft->AlignedSum(ir,adc);
            fft->AlignedSum(iw,signal);
      }
      if(geom->Plane(plane,tpc).SignalType() == geo::kCollection && 
	 rdvec[rd]->Channel() > colChan0 && 
	 rdvec[rd]->Channel() < colChan1) {
	fft->AlignedSum(cr,adc);
	fft->AlignedSum(cw,signal);
      }
      if(geom->Plane(plane,tpc).SignalType() == geo::kInduction) {
	if(*max_element(adc.begin(),adc.end()) > pedestal+threshold) 
	  fRawIndPeak->Fill(*max_element(adc.begin(),adc.end()));
	if(*max_element(signal.begin(),signal.end()) > pedestal+threshold) 
	  fCalIndPeak->Fill(*max_element(signal.begin(),signal.end()));
      }
      if(geom->Plane(plane,tpc).SignalType() == geo::kCollection) {
	if(*max_element(adc.begin(),adc.end()) > pedestal +threshold)
	  fRawColPeak->Fill(*max_element(adc.begin(),adc.end()));
	if(*max_element(signal.begin(),signal.end()) > pedestal+threshold) 
	  fCalColPeak->Fill(*max_element(signal.begin(),signal.end()));
      }
	
      int window = 8; 
      static unsigned int pulseHeight = 5;
      unsigned int tmin = 1; 
      unsigned int tmax = 1;
      int indMax = TMath::LocMax(signalSize,&adc[0]);
      double sigMin = 0.0;
      double sigMax = TMath::MaxElement(signalSize,&adc[0]);
      if(geom->Plane(plane,tpc).SignalType() == geo::kInduction && sigMax>=pulseHeight) {
	int indMin = TMath::LocMin(signalSize,&adc[0]);
	sigMin = TMath::MinElement(signalSize,&adc[0]);
	tmin = std::max(indMax-window,0);
	tmax = std::min(indMin+window,(int)signalSize-1);
	LOG_DEBUG("CalWireAna") << std::cout << "Induction channel, indMin,tmin,tmax " 
				   << rd<< " " << indMin<< " " << tmin << " " << tmax;
      } 
      else if (sigMax>=pulseHeight){
	tmin = std::max(indMax-window,0);
	tmax = std::min(indMax+window,(int)signalSize-1);
	LOG_DEBUG("CalWireAna") << "Collection channel, tmin,tmax "<< rd<< " " << tmin << " " << tmax;
      }
      
      fMin->Fill(sigMin);
      fMax->Fill(sigMax);
      fWindow->Fill(indMax, tmax - tmin);
      
      std::vector<double> winDiffs;
      int cnt=0;
      //	for(unsigned int t = tmin; t < tmax; ++t)
      static unsigned int tRawLead = 0;
      for(unsigned int t = 1; t < signalSize; ++t){	    
	fDiffsW->Fill(signal[t]-signal[t-1]);
	fWireSig->Fill(rd,signal[t]);
	
	if (t>=tmin && t<=tmax && tmax>=pulseHeight && (t+tRawLead)<signalSize){
	  cnt++;
	  winDiffs.push_back((adc[t+tRawLead]-signal[t])/adc[t+tRawLead]);
	  fDiffsRW->Fill(adc[t+tRawLead]-signal[t]);
	  fDiffsRWvsR->Fill(rd,adc[t+tRawLead]-signal[t]);
	  if (sigMax >= pulseHeight){
	    fDiffsRWgph->Fill(adc[t+tRawLead]-signal[t]);
	    fDiffsRWvsRgph->Fill(rd,adc[t+tRawLead]-signal[t]);
	  }
	}
      }               

      LOG_DEBUG("CalWireAna") << "on channel " << rdvec[rd]->Channel();
      // TMath::Mean with iterators doesn't work. EC,23-Sep-2010.
      double tmp = TMath::Mean(winDiffs.size(),&winDiffs[0]);
      double tmp2 = TMath::RMS(winDiffs.size(),&winDiffs[0]);
      double tmp3=0;
      for (int ii=0; ii<rdvec[rd]->Samples(); ii++) tmp3+=rdvec[rd]->ADC(ii);
      for(int i = 0; i < fft->FFTSize(); i++) {
	fIR->Fill(i,ir[i]);
	fIW->Fill(i,iw[i]);
	fCR->Fill(i,cr[i]);
	fCW->Fill(i,cw[i]);
      }
      fRD_WireMeanDiff2D->Fill(rd,tmp);
      fRD_WireRMSDiff2D->Fill(rd,tmp2);
      
    }//end loop over rawDigits

    return;
  }//end analyze method

}//end namespace
