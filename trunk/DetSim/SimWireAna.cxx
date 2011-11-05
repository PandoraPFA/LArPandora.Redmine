////////////////////////////////////////////////////////////////////////
// $Id: SimWireAna.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireAna class designed to make histograms
//
// brebel@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// LArSoft includes
#include "DetSim/SimWireAna.h"
#include "Geometry/geo.h"
#include "Simulation/sim.h"
#include "RawData/raw.h"

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
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

namespace detsim{

  //-------------------------------------------------
  SimWireAna::SimWireAna(fhicl::ParameterSet const& pset) 
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  SimWireAna::~SimWireAna()
  {
  }

  void SimWireAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fDetSimModuleLabel = p.get< std::string >("DetSimModuleLabel");

    return;
  }
  //-------------------------------------------------
  void SimWireAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fDiffs          = tfs->make<TH1F>("One timestamp diffs", ";#Delta ADC;",        40,   -19.5,  20.5);
    fCompressErr    = tfs->make<TH1F>("compressErr",         ";Raw-Compressed;",    1000, -495.5, 500.5);
    fCompressFactor = tfs->make<TH1F>("compressFactor",      ";Compression;",       500,     0.,    1.);

    fCompressErr2D  = tfs->make<TH2F>("compressErr2D",       ";Raw;Raw-Compressed", 100, -50., 50., 1000,  -495.5, 500.5);
    fRawVsCompress  = tfs->make<TH2F>("rawVsCompress",       ";Raw;Compressed",     100, -50., 50.,  100,   -50.,  50.);
  
    return;

  }

  //-------------------------------------------------
  void SimWireAna::analyze(const art::Event& evt)
  {

    // loop over the raw digits and get the adc vector for each, then compress it and uncompress it

    art::Handle< std::vector<raw::RawDigit> > rdHandle;
    evt.getByLabel(fDetSimModuleLabel,rdHandle);

    art::PtrVector<raw::RawDigit> rdvec;
    for(unsigned int i = 0; i < rdHandle->size(); ++i){
      art::Ptr<raw::RawDigit> r(rdHandle,i);
      rdvec.push_back(r);
    }

    /// loop over all the raw digits
    for(unsigned int rd = 0; rd < rdvec.size(); ++rd){

      std::vector<short> adc;
      std::vector<short> uncompressed(rdvec[rd]->Samples());
      for(unsigned int t = 1; t < rdvec[rd]->Samples(); ++t){
	fDiffs->Fill(rdvec[rd]->ADC(t) - rdvec[rd]->ADC(t-1));
	adc.push_back(rdvec[rd]->ADC(t-1));
      }
    
      //get the last one for the adc vector
      adc.push_back(rdvec[rd]->ADC(rdvec[rd]->Samples()-1));
    
      raw::Compress(adc, raw::kHuffman);
    
      fCompressFactor->Fill((1.*adc.size())/(1.*rdvec[rd]->Samples()));
    
      raw::Uncompress(adc, uncompressed, raw::kHuffman);
    
      if(uncompressed.size() != rdvec[rd]->Samples()){ 
	cet::exception("WrongSizeUncompress") 
	  << "uncompression does not produce same size vector as original: " 
	  << "original = " << rdvec[rd]->Samples() << " uncompress = " 
	  << uncompressed.size();
      }
    
      art::ServiceHandle<geo::Geometry> geo;
      unsigned int tpc           = 0;
      unsigned int plane         = 0;
      unsigned int wire          = 0;
      geo->ChannelToWire(rdvec[rd]->Channel(), tpc, plane, wire);
    
      for(unsigned int t = 0; t <  uncompressed.size(); ++t){
	//std::cout << t << " " << rdFE->ADC(t) << " " << uncompressed[t] << std::endl;
	if(uncompressed[t]-rdvec[rd]->ADC(t) > 1) std::cout << "problem with event " //<<evt
							    << " time " << t << " ADC " << rdvec[rd]->ADC(t) 
							    << " uncompress " << uncompressed[t] 
							    << " tpc " << tpc
							    << " plane " << plane 
							    << " wire " << wire << std::endl;
	
	fCompressErr->Fill(uncompressed[t]-rdvec[rd]->ADC(t));
	fCompressErr2D->Fill(rdvec[rd]->ADC(t), uncompressed[t]-rdvec[rd]->ADC(t));
	fRawVsCompress->Fill(rdvec[rd]->ADC(t), uncompressed[t]);
      }
    }//end loop over digits

    return;
  }//end analyze method

}//end namespace
