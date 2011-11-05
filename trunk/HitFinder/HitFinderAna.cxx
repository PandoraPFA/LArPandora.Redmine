////////////////////////////////////////////////////////////////////////
// $Id: HitFinderAna.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// HitFinderAna class designed to make histograms
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
#include "HitFinder/HitFinderAna.h"
#include "Geometry/geo.h"
#include "MCCheater/BackTracker.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"

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
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace hit{

  //-------------------------------------------------
  HitFinderAna::HitFinderAna(fhicl::ParameterSet const& pset) 
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  HitFinderAna::~HitFinderAna()
  {
  }

  void HitFinderAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fFFTHitFinderModuleLabel = p.get< std::string >("HitsModuleLabel");
    fLArG4ModuleLabel        = p.get< std::string >("LArGeantModuleLabel");
    return;
  }
  //-------------------------------------------------
  void HitFinderAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    fNp0 = 9000;
    fNp1 = 9000;
    fNp2 = 9000;

    fHTree = tfs->make<TTree>("HTree","HTree");
    fTimep0 = new Float_t[fNp0];
    fTimep1 = new Float_t[fNp1];
    fTimep2 = new Float_t[fNp2];
    fWirep0 = new Int_t[fNp0];
    fWirep1 = new Int_t[fNp1];
    fWirep2 = new Int_t[fNp2];
    fChgp0 = new Float_t[fNp0];
    fChgp1 = new Float_t[fNp1];
    fChgp2 = new Float_t[fNp2];
    fXYZp0 = new Float_t[fNp0*3];
    fXYZp1 = new Float_t[fNp1*3];
    fXYZp2 = new Float_t[fNp2*3];

    fMCPdg0 = new Int_t[fNp0];
    fMCPdg1 = new Int_t[fNp1];
    fMCPdg2 = new Int_t[fNp2];
    fMCTId0 = new Int_t[fNp0];
    fMCTId1 = new Int_t[fNp1];
    fMCTId2 = new Int_t[fNp2];
    fMCE0 = new Float_t[fNp0];
    fMCE1 = new Float_t[fNp1];
    fMCE2 = new Float_t[fNp2];

    fHTree->Branch("HEvt", &fEvt, "HEvt/I");
    fHTree->Branch("HRun", &fRun, "HRun/I");
    fHTree->Branch("HNp0", &fNp0, "HNp0/I");
    fHTree->Branch("HNp1", &fNp1, "HNp1/I");
    fHTree->Branch("HNp2", &fNp2, "HNp2/I");
    fHTree->Branch("HN3p0", &fN3p0, "HN3p0/I");
    fHTree->Branch("HN3p1", &fN3p1, "HN3p1/I");
    fHTree->Branch("HN3p2", &fN3p2, "HN3p2/I");
    fHTree->Branch("Htp0", fTimep0, "Htp0[HNp0]/F");
    fHTree->Branch("Htp1", fTimep1, "Htp1[HNp1]/F");
    fHTree->Branch("Htp2", fTimep2, "Htp2[HNp2]/F");
    fHTree->Branch("Hwp0", fWirep0, "Hwp0[HNp0]/I");
    fHTree->Branch("Hwp1", fWirep1, "Hwp1[HNp1]/I");
    fHTree->Branch("Hwp2", fWirep2, "Hwp2[HNp2]/I");
    fHTree->Branch("Hchgp0", fChgp0, "Hchgp0[HNp0]/F");
    fHTree->Branch("Hchgp1", fChgp1, "Hchgp1[HNp1]/F");
    fHTree->Branch("Hchgp2", fChgp2, "Hchgp2[HNp2]/F");
    fHTree->Branch("HMCXYZp0", fXYZp0, "HMCXYZp0[HN3p0]/F");
    fHTree->Branch("HMCXYZp1", fXYZp1, "HMCXYZp1[HN3p1]/F");
    fHTree->Branch("HMCXYZp2", fXYZp2, "HMCXYZp2[HN3p2]/F");
    fHTree->Branch("HMCPdgp0", fMCPdg0, "HMCPdgp0[HNp0]/I");
    fHTree->Branch("HMCPdgp1", fMCPdg1, "HMCPdgp1[HNp1]/I");
    fHTree->Branch("HMCPdgp2", fMCPdg2, "HMCPdgp2[HNp2]/I");
    fHTree->Branch("HMCTIdp0", fMCTId0, "HMCTIdp0[HNp0]/I");
    fHTree->Branch("HMCTIdp1", fMCTId1, "HMCTIdp1[HNp1]/I");
    fHTree->Branch("HMCTIdp2", fMCTId2, "HMCTIdp2[HNp2]/I");
    fHTree->Branch("HMCEp0", fMCE0, "HMCEp0[HNp0]/F");
    fHTree->Branch("HMCEp1", fMCE1, "HMCEp1[HNp1]/F");
    fHTree->Branch("HMCEp2", fMCE2, "HMCEp2[HNp2]/F");

  
    return;

  }

  //-------------------------------------------------
  void HitFinderAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()){
      throw cet::exception("HitFinderAna: ") << "Not for use on Data yet... " << "\n";
    }
    
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fFFTHitFinderModuleLabel,hitHandle);

    sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArG4ModuleLabel);
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fLArG4ModuleLabel, sccol);

    std::cout << _particleList << std::endl;

    //    art::PtrVector<recob::Hit> hits;
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitHandle);
    
    art::ServiceHandle<geo::Geometry> geom;  
  
    unsigned int p(0),w(0), t(0), channel(0);
    for(unsigned int tpc = 0; tpc < geom->NTPC(); ++tpc){
      //      for(unsigned int plane=0; plane<geom->Nplanes(tpc); plane++){
      //	for(unsigned int i = 0; i< hitcol->size(); ++i){
      fNp0=0;       fN3p0=0;
      fNp1=0;       fN3p1=0;
      fNp2=0;       fN3p2=0;

      //now make a vector where each channel in the detector is an entry
      std::vector<const sim::SimChannel*> scs(geom->Nchannels(),0);
      for(size_t i = 0; i < sccol.size(); ++i) scs[sccol[i]->Channel()] = sccol[i];

      std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
      while(itr != hits.end()) {

	//art::Ptr<recob::Hit> hit(hitHandle, i);
	channel=(*itr)->Wire()->RawDigit()->Channel();
	t=0;p=0;w=0;
	geom->ChannelToWire(channel,t,p,w);

	fRun = evt.run();
	fEvt = evt.id().event();

	/*
	std::cout << "HitFinderAna: channel # is " << channel << std::endl;
	std::cout << "HitFinderAna: Hit itself is " << (*itr) << std::endl;
	std::cout << "HitFinderAna: SimIDE itself is " << scs[channel] << std::endl;
	*/
	if (!scs[channel]) {itr++;continue;}

	std::vector<cheat::TrackIDE> trackides = cheat::BackTracker::HitToTrackID(*(scs[(*itr)->Channel()]), *itr);
	std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
	std::vector<double> xyz = cheat::BackTracker::HitToXYZ(*(scs[(*itr)->Channel()]),*itr);


	if (p==0 && fNp0<9000) 
	  {
	    fTimep0[fNp0] = (*itr)->PeakTime();
	    fWirep0[fNp0] = w;
	    fChgp0[fNp0] = (*itr)->Charge();

	    for (unsigned int kk=0;kk<3;kk++)
	      {
		fXYZp0[fNp0*3+kk] = xyz[kk];
	      }
	    

	    while( idesitr != trackides.end() )
	      {
		fMCTId0[fNp0] = (*idesitr).trackID;
		if (_particleList.find((*idesitr).trackID) != _particleList.end()) 
		  {
		    const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
		    fMCPdg0[fNp0] = particle->PdgCode();
		    fMCE0[fNp0] = particle->E();
		  }
		idesitr++;
	      }

	    fNp0++;
	  }

	else if (p==1 && fNp1<9000) 
	  {
	    fTimep1[fNp1] = (*itr)->PeakTime();
	    fWirep1[fNp1] = w;
	    fChgp1[fNp1] = (*itr)->Charge();

	    for (unsigned int kk=0;kk<3;kk++)
	      {
		fXYZp1[fNp1*3+kk] = xyz[kk];
	      }

	    while( idesitr != trackides.end() )
	      {
		fMCTId1[fNp1] = (*idesitr).trackID;
		if (_particleList.find((*idesitr).trackID) != _particleList.end())
		  {
		    const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
		    fMCPdg1[fNp1] = particle->PdgCode();
		    fMCE1[fNp1] = particle->E();
		  }
		idesitr++;
	      }
	    fNp1++;
	  }

	else if (p==2  && fNp2<9000) 
	  {
	    fTimep2[fNp2] = (*itr)->PeakTime();
	    fWirep2[fNp2] = w;
	    fChgp2[fNp2] = (*itr)->Charge();
	    
	    for (unsigned int kk=0;kk<3;kk++)
	      {
		fXYZp2[fNp2*3+kk] = xyz[kk];
	      }
	    
	    while( idesitr != trackides.end())
	      {
		fMCTId2[fNp2] = (*idesitr).trackID;
		if (_particleList.find((*idesitr).trackID) != _particleList.end() ) 
		  {
		    const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
		    fMCPdg2[fNp2] = particle->PdgCode();
		    fMCE2[fNp2] = particle->E();
		  }
		idesitr++;
	      }
	    fNp2++;
	  }

	fN3p0 = 3* fNp0;
	fN3p1 = 3* fNp1;
	fN3p2 = 3* fNp2;

	fHTree->Fill();
	itr++;
      } // loop on Hits
      //      }
    } //  loop on NTPCs

    return;
  }//end analyze method

}//end namespace
