////////////////////////////////////////////////////////////////////////
//
// Kinematics class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to pull information from the event record and put it together coherently
////////////////////////////////////////////////////////////////////////


#include <iostream>


// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

#include "T962/Kinematics/Kinematics.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"
#include <TTree.h>

#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/LArVoxelCalculator.h"
#include "Simulation/LArVoxelData.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "T962/T962_Objects/MINOS.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

namespace kin {
//-----------------------------------------------------------------------------
kin::Kinematics::Kinematics(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel            (pset.get< std::string >("DBScanModuleLabel")),
  fLArG4ModuleLabel             (pset.get< std::string >("LArG4ModuleLabel")),
  fGenieGenModuleLabel             (pset.get< std::string >("GenieGenModuleLabel")),
  fScanModuleLabel              (pset.get< std::string > ("ScanModuleLabel")),
  fMINOSModuleLabel              (pset.get< std::string > ("MINOSModuleLabel")),
  fLarTracks_label               (pset.get< std::string >("LArTracksModuleLabel")),
  fm_event(0),
  fm_pdgcode(0),
  fm_leptheta_reco(0.),
  fm_leptheta_true(0.),
  fm_lepphi_reco(0.),
  fm_lepphi_true(0.),
  fm_init_x_reco(-99.),
  fm_init_y_reco(-99.),
  fm_init_z_reco(-99.),
  fm_init_x_true(-99.),
  fm_init_y_true(-99.),
  fm_init_z_true(-99.),
  fm_tpcexit_x_reco(-99.),
  fm_tpcexit_y_reco(-99.),
  fm_tpcexit_z_reco(-99.),
  fm_tpcexit_x_true(-99.),
  fm_tpcexit_y_true(-99.),
  fm_tpcexit_z_true(-99.)  
{

}

//-----------------------------------------------------------------------------
Kinematics::~Kinematics()
{
}

//-----------------------------------------------------------------------------


void Kinematics::beginJob()
{
  // get access to the TFile service  
  art::ServiceHandle<art::TFileService> tfs;  
  ftree= tfs->make<TTree>("KinTree","KinTree");
  ftree->Branch("event", &fm_event, "event/I");
  ftree->Branch("pdgcode", &fm_pdgcode, "pdgcode/I");
  ftree->Branch("leptheta_reco", &fm_leptheta_reco, "leptheta_reco/F");
  ftree->Branch("leptheta_true", &fm_leptheta_true, "leptheta_true/F");
  ftree->Branch("lepphi_reco", &fm_lepphi_reco, "lepphi_reco/F");
  ftree->Branch("lepphi_true", &fm_lepphi_true, "lepphi_true/F");
  ftree->Branch("init_x_reco", &fm_init_x_reco, "init_x_reco/F");
  ftree->Branch("init_y_reco", &fm_init_y_reco, "init_y_reco/F");
  ftree->Branch("init_z_reco", &fm_init_z_reco, "init_z_reco/F");
  ftree->Branch("init_x_true", &fm_init_x_true, "init_x_true/F");
  ftree->Branch("init_y_true", &fm_init_y_true, "init_y_true/F");
  ftree->Branch("init_z_true", &fm_init_z_true, "init_z_true/F");
  ftree->Branch("tpcexit_x_reco", &fm_tpcexit_x_reco, "tpcexit_x_reco/F");
  ftree->Branch("tpcexit_y_reco", &fm_tpcexit_y_reco, "tpcexit_y_reco/F");
  ftree->Branch("tpcexit_z_reco", &fm_tpcexit_z_reco, "tpcexit_z_reco/F");
  ftree->Branch("tpcexit_x_true", &fm_tpcexit_x_true, "tpcexit_x_true/F");
  ftree->Branch("tpcexit_y_true", &fm_tpcexit_y_true, "tpcexit_y_true/F");
  ftree->Branch("tpcexit_z_true", &fm_tpcexit_z_true, "tpcexit_z_true/F");
}

//-----------------------------------------------------------------------------
void Kinematics::analyze(const art::Event& evt) 
{
  double vertex [3] = { 0, 0, 0 };
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larp;
  double electronlifetime=larp->ElectronLifetime();
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);

  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  
  art::Handle< std::vector<recob::Track> > LarTrackHandle;
  evt.getByLabel(fLarTracks_label,LarTrackHandle);
  double larStart[3]={0,0,0};
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
 
  double larEnd[3];
     //just take the longest 3D track for now    
    //      for(unsigned int i=0; i<LarTrackHandle->size();++i){
    
    // std::sort(LarTrackHandle->begin(),LarTrackHandle->end());
    int size=0;
     if(LarTrackHandle->size()>0)
     for(unsigned int i=0; i<LarTrackHandle->size();++i){
      art::Ptr<recob::Track> lartrack(LarTrackHandle,i);
      
       std::cout << " T962 " << *lartrack << std::endl;
       std::cout << "T962 Track #" << lartrack->ID() <<std::endl;

       if(lartrack->SpacePoints().size()>size)
       {
       lartrack->Direction(larStart,larEnd);
       lartrack->Extent(trackStart,trackEnd);  
       }
       size=lartrack->SpacePoints().size();
       std::cout<<"size "<<size<<std::endl;
      }
      
   if(LarTrackHandle->size()==0)
   std::cout << " *************************Failed on event " <<evt.id().event()<< std::endl;
      
   fm_event=evt.id().event(); 


   art::PtrVector<simb::MCTruth> mclist;
   for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 


    for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
    simb::MCParticle neut(mc->GetParticle(i));

    if(mclist.size()>1)
    break;

    fm_init_x_true=neut.Vx();
    fm_init_y_true=neut.Vy();
    fm_init_z_true=neut.Vz();
    
    
    
    fm_leptheta_true=TMath::ACos(neut.Pz()/sqrt(pow(neut.Px(),2)+pow(neut.Py(),2)+pow(neut.Pz(),2)));
    
    fm_lepphi_true=(TMath::Pi()+TMath::ATan2(neut.Py(),neut.Px()));

    fm_lepphi_reco=TMath::Pi()+TMath::ATan2(larStart[1],larStart[0]);
    fm_leptheta_reco=TMath::ACos(larStart[2]);
    if(trackStart.size())
    {
    fm_init_x_reco=trackStart[0];
    fm_init_y_reco=trackStart[1];
    fm_init_z_reco=trackStart[2];
    }
    else
    {
    fm_init_x_reco=-99.;
    fm_init_y_reco=-99.;
    fm_init_z_reco=-99.;    
    }
    
    if(trackEnd.size())
    {
    fm_tpcexit_x_reco=trackEnd[0];
    fm_tpcexit_y_reco=trackEnd[1];
    fm_tpcexit_z_reco=trackEnd[2];
    }
    else
    {
    fm_tpcexit_x_reco=-99.;
    fm_tpcexit_y_reco=-99.;
    fm_tpcexit_z_reco=-99.;    
    }
   }

    // get the particles from the event handle
    art::Handle< std::vector<sim::Particle> > parHandle;
    evt.getByLabel(fLArG4ModuleLabel, parHandle);

    art::PtrVector<sim::Particle> pvec;
    for(unsigned int i = 0; i < parHandle->size(); ++i){
      art::Ptr<sim::Particle> p(parHandle, i);      
      pvec.push_back(p);
    }
    
    for(unsigned int i = 0; i < pvec.size(); ++i)
    {
     if(pvec[i]->Process()!="primary")
     continue; 
    simb::MCTrajectory trajectory = pvec[i]->Trajectory();
    fm_pdgcode=pvec[i]->PdgCode();       
    int numberofpoints= pvec[i]->NumberTrajectoryPoints();   

     for(int j=1; j<numberofpoints; j++)
     {
      TLorentzVector prevposition = trajectory.Position(j-1);
      TLorentzVector position = trajectory.Position(j);
      TLorentzVector prevmomentum = trajectory.Momentum(j-1);
      TLorentzVector momentum = trajectory.Momentum(j);

       if(       (((prevposition.X()<0.&&position.X()>=0.)||(prevposition.X()>0.&&position.X()<=0.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.X()<47.&&position.X()>=47.)||(prevposition.X()>47.&&position.X()<=47.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.Y()<20.&&position.Y()>=20.)||(prevposition.Y()>20.&&position.Y()<=20.))&&(prevposition.X()<47.&&prevposition.X()>=0.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.Y()<-20.&&position.Y()>=-20.)||(prevposition.Y()>-20.&&position.Y()<=-20.))&&(prevposition.X()<47.&&prevposition.X()>=0.&&prevposition.Z()>0.&&prevposition.Z()<=90.))
       ||       (((prevposition.Z()<0.&&position.Z()>=0.)||(prevposition.Z()>0.&&position.Z()<=0.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.X()>0.&&prevposition.X()<=47.))
       ||       (((prevposition.Z()<90.&&position.Z()>=90.)||(prevposition.Z()>90.&&position.Z()<=90.))&&(prevposition.Y()<20.&&prevposition.Y()>=-20.&&prevposition.X()>0.&&prevposition.X()<=47.))
       )
       {    
        fm_tpcexit_x_true=position.X();
        fm_tpcexit_y_true=position.Y();
        fm_tpcexit_z_true=position.Z();
       }
       
       std::cout<<fm_pdgcode<<" at z position "<<position.Z()<<" lost "<<((prevmomentum.E()-momentum.E())*1000.)<<"MeV in "<<(position.Z()-prevposition.Z())<<"cm which is "<<((prevmomentum.E()-momentum.E())*1000.)/(position.Z()-prevposition.Z())<< "MeV/cm VolumeName is " << geom->VolumeName(position.Vect())<<" "<<geom->MaterialName(position.Vect()) << std::endl;  
       
     }
    }

ftree->Fill(); 

//hand scan info
//    art::PtrVector<t962::ScanInfo> scanIn;
//    scanIn.clear();
//    art::Handle< std::vector<t962::ScanInfo> > scanHandle;
//    evt.getByLabel(fScanModuleLabel,scanHandle);
//  
//   for(unsigned int i = 0; i < scanHandle->size(); ++i){     
//       art::Ptr<t962::ScanInfo> scaninfo(scanHandle, i);
//        scanIn.push_back(scaninfo);     
//      }
//  
//      for(unsigned int i = 0; i < scanIn.size(); ++i){     
// //       vertexcolwire=scanIn[i]->Get_VertColWire();
// //       vertexcoltime=scanIn[i]->Get_VertColTime();   
//      }
     
//      std::cout<<"vertex: "<<vertexcolwire<<" "<<vertexcoltime<<std::endl;
//      findcol->Fill(vertexcolwire,vertexcoltime);
//      
     
 //minos info
//    art::PtrVector<t962::MINOS> minosIn;
//    minosIn.clear();
//    art::Handle< std::vector<t962::MINOS> > minosHandle;
//    evt.getByLabel(fMINOSModuleLabel,minosHandle);
//  
//   for(unsigned int i = 0; i < minosHandle->size(); ++i){     
//       art::Ptr<t962::MINOS> minosinfo(minosHandle, i);
//        minosIn.push_back(minosinfo);     
//      }
//  
//      for(unsigned int i = 0; i < minosIn.size(); ++i){     
//       std::cout<<"trkd cosy: "<<minosIn[i]->ftrkdcosy<<std::endl;
//    
//      }
//      
    
    
//take into account automated vertex finding as well (not yet implemented)  
//    art::PtrVector<recob::Vertex> vertIn;
//  
//    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
//      {
//        art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
//        vertIn.push_back(vertex);
//      }
//    
//     for(unsigned int i = 0; i < vertIn.size(); ++i){     
//       std::cout<<vertIn[i]->WireNum()<<" "<<vertIn[i]->DriftTime()<<std::endl;
//      }  

//   filter::ChannelFilter chanFilt;  
//   art::PtrVector<recob::Hit> cHits;
//   art::PtrVector<recob::Hit> hit;
//    
//   art::PtrVector<recob::Cluster> clusIn;
//   for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
//     {
//       art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
//       clusIn.push_back(cluster);
//     }

 }

}//namespace
