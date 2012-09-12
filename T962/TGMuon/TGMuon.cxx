////////////////////////////////////////////////////////////////////////
/// \file  TGMuon.cxx
/// \brief Generator for through-going muons
/// Module designed to produce through-going muons based on MINOS data for ArgoNeuT (run 3, neutrino-mode right now)
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <dirent.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>

// nutools includes
#include "SimulationBase/simbase.h"

// lar includes
#include "T962/TGMuon/TGMuon.h"
#include "Geometry/geo.h"
#include "SummaryData/summary.h"

#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TRandom3.h"

namespace t962{

  //____________________________________________________________________________
  TGMuon::TGMuon(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
    produces< sumdata::POTSummary, art::InSubRun >();
    fRand = new TRandom3();
    fRand->SetSeed();
    totalpot=0.;    
  }

  //____________________________________________________________________________
  TGMuon::~TGMuon()
  {
  if(fRand) delete fRand;
  }

  //____________________________________________________________________________

  //____________________________________________________________________________
  void TGMuon::beginRun(art::Run& run)
  {
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    geo::DetId_t detid = geo->DetId();
    std::auto_ptr<sumdata::RunData> runcol(new sumdata::RunData(detid));
    run.put(runcol);
    return;
  }

  void TGMuon::endSubRun(art::SubRun& sr)
  {
    std::auto_ptr<sumdata::POTSummary> p(new sumdata::POTSummary());    
    p->totpot = totalpot;
    p->totgoodpot = totalpot;    
    sr.put(p);    
    return;
  }

  //____________________________________________________________________________
  void TGMuon::produce(art::Event& evt)
  {
    ///auto_ptr allows ownership to be transferred to the art::Event after the put statement
    std::auto_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    Sample(truth);
    truthcol->push_back(truth);
    evt.put(truthcol);

    return;
  }

  //____________________________________________________________________________
  void TGMuon::Sample(simb::MCTruth &mct) 
  {
      //loop through MINOS root files	  
      std::string file = "/argoneut/data/outstage/spitz7/allminos_neutrinodata/allminos_neutrinodata.root";        
      int snarl=0;
      int run=0;
      int subRun=0;
      int flag=0;
      double r=0;
      double x_offset=117.4; 
      double y_offset=19.3;
      double z_offset=150.; //force particle to be well upstream of ArgoNeuT
  
            TFile *f = new TFile(file.c_str());
            TTree *minitree = (TTree*)f->Get("minitree");
            minitree->SetBranchAddress("run",&frun);  
            minitree->SetBranchAddress("subRun",&fsubRun);            
            minitree->SetBranchAddress("trkIndex",&ftrkIndex);
            minitree->SetBranchAddress("ntrkstp",&fntrkstp);
            minitree->SetBranchAddress("snarl",&fsnarl);
            minitree->SetBranchAddress("trtgtd",&ftrtgtd);
            minitree->SetBranchAddress("trkVtxX",&ftrkVtxX);
            minitree->SetBranchAddress("trkVtxY",&ftrkVtxY);
            minitree->SetBranchAddress("trkVtxZ",&ftrkVtxZ);
            minitree->SetBranchAddress("trkdcosx",&ftrkdcosx);
            minitree->SetBranchAddress("trkdcosy",&ftrkdcosy);
            minitree->SetBranchAddress("trkdcosz",&ftrkdcosz);
            minitree->SetBranchAddress("trkmom",&ftrkmom);
            minitree->SetBranchAddress("trkqp",&ftrkqp);            
            Long64_t nentries = minitree->GetEntries();
            //Long64_t nbytes = 0;           
      
      int flag2=0;
      while(flag2==0)    
      {      
            r=(int)fRand->Uniform(1,nentries-1);
            minitree->GetEntry(r);
            run=frun;
            subRun=fsubRun;
            snarl=fsnarl;
            //count POT even if the particle is not passed to Geant
            totalpot+=ftrtgtd*pow(10,12);
          
      for(int i=r-ftrkIndex;i<nentries;i++)
      {
          minitree->GetEntry(i);
           
          if(frun==run && fsubRun==subRun && fsnarl==snarl)
          flag=1;           
           
    	  if(flag==1 && fsnarl!=snarl)
     	  break;
        
     	 TVector3 x;
		 x[0]=100.*ftrkVtxX - x_offset;
    	 x[1]=100.*ftrkVtxY + y_offset;
    	 x[2]=100.*ftrkVtxZ - z_offset;        
        //don't bother passing the particle to geant if it's way outside the argoneut detector window   
      	if(100.0*ftrkVtxZ > 20. || fabs(x[0]+23.5)>100. || fabs(x[1])>100.)
      	continue;

      	// Choose momentum
      	double p = 0.0;
      	double m = 0.0;
      	p=ftrkmom;
      
      	if(!p)
      	continue;
      
      	int pdg=0;
      	
      	if(ftrkqp>0)
      	pdg=-13;
      	else
     	pdg=13;
      
      	TParticlePDG* pdgp = TDatabasePDG::Instance()->GetParticle(pdg);
      	if (pdgp) m = pdgp->Mass();

      	TLorentzVector pos(x[0], x[1], x[2], 0.0);
      	TLorentzVector pvec(ftrkdcosx*p,ftrkdcosy*p,ftrkdcosz*p,sqrt(p*p+m*m));
    
      	int trackid = -1*(i+1); // set track id to -i as these are all primary particles and have id <= 0
      	std::string primary("primary");

      	simb::MCParticle part(trackid, pdg, primary);
      	part.AddTrajectoryPoint(pos, pvec);       
      	mct.Add(part);
	flag2=1;      
      }//for loop
      }//while flag2 loop
       minitree->Delete();
       f->Close();
       f->Delete();

   return;
  }

}//end namespace evgen
////////////////////////////////////////////////////////////////////////
