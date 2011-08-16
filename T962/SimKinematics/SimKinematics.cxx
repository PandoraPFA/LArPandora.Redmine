////////////////////////////////////////////////////////////////////////
//
// VertexActivity class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to characterize the vertex activity associated with a neutrino event
////////////////////////////////////////////////////////////////////////


#include <iostream>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/OrphanHandle.h"


#include "T962/SimKinematics/SimKinematics.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "T962/T962_Objects/MINOS.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"


//-----------------------------------------------------------------------------
simkin::SimKinematics::SimKinematics(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel            (pset.get< std::string >("DBScanModuleLabel")),
  fLArG4ModuleLabel             (pset.get< std::string >("LArG4ModuleLabel")),
  fGenieGenModuleLabel             (pset.get< std::string >("GenieGenModuleLabel")),
  fScanModuleLabel              (pset.get< std::string > ("ScanModuleLabel")),
  fMINOSModuleLabel              (pset.get< std::string > ("MINOSModuleLabel")),
  fm_run(0), 
  fm_event(0), 
  fm_CCNC(0),
  fm_PDG(0),
  fm_mode(0),
  fm_hitnuc(0),
  fm_leppx(0.),
  fm_leppy(0.),
  fm_leppz(0.),
  fm_lepE(0.),
  fm_W(0.),
  fm_qsqr(0.),
  fm_vertx(0.),
  fm_verty(0.),
  fm_vertz(0.),
  fm_lepphi(0.),
  fm_leptheta(0.)
{
  //produces< std::vector<recob::Vertex> >();
}

//-----------------------------------------------------------------------------
simkin::SimKinematics::~SimKinematics()
{
}

//-----------------------------------------------------------------------------



void simkin::SimKinematics::beginJob()
{
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  
  ftree= tfs->make<TTree>("SimKinTree","SimKinTree");
  ftree->Branch("run", &fm_run, "run/I");
  ftree->Branch("event", &fm_event, "event/I");
  ftree->Branch("CCNC", &fm_CCNC, "CCNC/I");
  ftree->Branch("PDG", &fm_PDG, "PDG/I");
  ftree->Branch("mode", &fm_mode, "mode/I");
  ftree->Branch("hitnuc", &fm_hitnuc, "hitnuc/I");

  ftree->Branch("leppx", &fm_leppx, "leppx/F");
  ftree->Branch("leppy", &fm_leppy, "leppy/F");
  ftree->Branch("leppz", &fm_leppz, "leppz/F");
  ftree->Branch("lepphi", &fm_lepphi, "lepphi/F");
  ftree->Branch("leptheta", &fm_leptheta, "leptheta/F");
  ftree->Branch("lepE", &fm_lepE, "lepE/F");
  ftree->Branch("W", &fm_W, "W/F");
  ftree->Branch("qsqr", &fm_qsqr, "qsqr/F");
  ftree->Branch("vertx", &fm_vertx, "vertx/F");
  ftree->Branch("verty", &fm_verty, "verty/F");
  ftree->Branch("vertz", &fm_vertz, "vertz/F");
}

//-----------------------------------------------------------------------------
void simkin::SimKinematics::produce(art::Event& evt)
{
  double vertex [3] = { 0, 0, 0 };
  art::ServiceHandle<geo::Geometry> geom;
  
  std::cout<<"half width geom->DetHalfWidth() "<<geom->DetHalfWidth()<<std::endl;
  std::cout<<"half height geom->DetHalfHeight() "<<geom->DetHalfHeight()<<std::endl;
  std::cout<<"length geom->DetLength() "<<geom->DetLength()<<std::endl;
  
  art::ServiceHandle<util::LArProperties> larp;
  double electronlifetime=larp->ElectronLifetime();
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);

TVector3 pxpypz;
//      
    
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);

  art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
   

   
   
  fm_run=evt.id().run();
  fm_event=evt.id().event(); 
   
  //  neutrinos  
//      for( unsigned int i = 0; i < mclist.size(); ++i ){

     art::Ptr<simb::MCTruth> mc(mclist[0]);

	simb::MCNeutrino neut(mc->GetNeutrino());



    fm_vertx =neut.Nu().Vx();
    fm_verty =neut.Nu().Vy();
    fm_vertz =neut.Nu().Vz();
    fm_CCNC=neut.CCNC();
    fm_PDG=neut.Lepton().PdgCode();
    fm_mode=neut.Mode();
    fm_hitnuc=neut.HitNuc();        
    fm_leppx=neut.Lepton().Px();
    fm_leppy=neut.Lepton().Py();
    fm_leppz=neut.Lepton().Pz();
    fm_lepE=neut.Lepton().E();
    fm_nuE=neut.Nu().E();
    fm_W=neut.W();
    fm_qsqr=neut.QSqr();
    fm_lepphi=(TMath::Pi()+TMath::ATan2(neut.Lepton().Py(),neut.Lepton().Px()));
    
    fm_leptheta=TMath::ACos(neut.Lepton().Pz()/sqrt(pow(neut.Lepton().Px(),2)+pow(neut.Lepton().Py(),2)+pow(neut.Lepton().Pz(),2)));
    
    
    ftree->Fill();     
  
// std::cout<<"Energy "<<neut.Nu().E()<<" "<<neut.CCNC()<<" "<<neut.Lepton().E()<<" "<<(180./3.141)*TMath::ACos(neut.Lepton().Pz()/sqrt(pow(neut.Lepton().Px(),2)+pow(neut.Lepton().Py(),2)+pow(neut.Lepton().Pz(),2)))<<" "<<neut.Mode()<<" "<<neut.HitNuc()<<std::endl;
// 
// 
// if(neut.CCNC()==0&&neut.Lepton().PdgCode()==13)
// {
// 
// lepangle->Fill((180./3.141)*TMath::ACos(neut.Lepton().Pz()/sqrt(pow(neut.Lepton().Px(),2)+pow(neut.Lepton().Py(),2)+pow(neut.Lepton().Pz(),2))));
// lepenergy->Fill(neut.Lepton().E());
// }


// }

 }



