////////////////////////////////////////////////////////////////////////
//
// CCQE MC analyzer
//
// \author kinga.partyka@yale.edu
// 
// 
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TMath.h"

#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "T962/CCQEanalysis/CCQEanalysis.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"
#include "Simulation/LArVoxelCalculator.h"
#include "Simulation/LArVoxelData.h"
#include "Simulation/LArVoxelID.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"


 
//-------------------------------------------------
t962::CCQEanalysis::CCQEanalysis(fhicl::ParameterSet const& pset) : 
 
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fClusterFinderModuleLabel (pset.get< std::string >("ClusterFinderModuleLabel"))
{



}

//-------------------------------------------------
t962::CCQEanalysis::~CCQEanalysis()
{

}

void t962::CCQEanalysis::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
Mu_theta=tfs->make<TH1F>("Mu_theta","Muon theta angle", 360,0 ,360);
Mu_phi=tfs->make<TH1F>("Mu_phi","Muon phi angle", 360,0 ,360);
mu_charge=tfs->make<TH1F>("mu_charge","Muon charge", 4,-2 ,2);
No_protons_in_event=tfs->make<TH1F>("No_protons_in_event","No of protons in each CCQE event", 10,0 ,10);
No_particles_in_event=tfs->make<TH1F>("No_particles_in_event","No of particles in each CCQE event (excluding neutrons), track>1 cm", 20,0 ,20);
P_containment=tfs->make<TH1F>("P_containment","Proton contained=1, else 0", 3,-1 ,2);
proton_track_length=tfs->make<TH1F>("proton_track_length","Length of proton track", 300,0 ,150);

Vertex_x=tfs->make<TH1F>("Vertex_x","vertex X coordinate value", 60,0 ,60);
Vertex_y=tfs->make<TH1F>("Vertex_y","vertex Y coordinate value", 60,-30 ,30);
Vertex_z=tfs->make<TH1F>("Vertex_z","vertex Z coordinate value", 100,0 ,100);

  
}

void t962::CCQEanalysis::analyze(const art::Event& evt)
{
  
  std::cout << "run    : " << evt.run() << std::endl;
  std::cout << "event  : " << evt.id().event() << std::endl;
  //----------------------------------------------------------------

 
  if (evt.isRealData()) 
    {
      std::cout<<"**** Don't call this module if you're not MC. "<<std::endl;
      return;
    }

  
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  
  art::Handle< std::vector<sim::LArVoxelData> > vxlistHandle;
  evt.getByLabel(fLArG4ModuleLabel,vxlistHandle);
  
  sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArG4ModuleLabel);
  
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
 
   art::Handle< std::vector<sim::Particle> > geant_list;
    evt.getByLabel (fLArG4ModuleLabel,geant_list);
    
  
  
  
   //................................................................
double vertex [3] = { 0, 0, 0 };
int have_p=0;
int have_pion=0;
double E_p, E_pion;
int event_has_pi0=0;
int event_has_pi_plus=0;

double MC_Total_Eng=0;

   art::PtrVector<simb::MCTruth> mclist;
   for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
    
     art::PtrVector<sim::Particle> geant_part;
   for (unsigned int ii = 0; ii <  geant_list->size(); ++ii)
    {
      art::Ptr<sim::Particle> p(geant_list,ii);
      geant_part.push_back(p);
    } 
    
    
    
    
    
  //..........................................................
   
    //First determine what kind of event we are dealing with, and select only the ones you want:
 //..........................................................
 int no_protons=0;
 
 
 for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
    simb::MCParticle neut(mc->GetParticle(i));
    
    vertex[0] =neut.Vx();
    vertex[1] =neut.Vy();
    vertex[2] =neut.Vz();
    std::cout<<"neut.Vx()= "<<neut.Vx()<<" ,y= "<<neut.Vy()<<" ,z= "<<neut.Vz()<<std::endl;
    
    Vertex_x->Fill(neut.Vx());
    Vertex_y->Fill(neut.Vy());
    Vertex_z->Fill(neut.Vz());
    std::cout<<"List from Genie: "<<std::endl;
     for(int j = 0; j < mc->NParticles(); ++j){
    simb::MCParticle part(mc->GetParticle(j));
    
    std::cout<<"pdg= "<<part.PdgCode()<<" ,Process="<<part.Process()<<" StatusCode= "<<part.StatusCode()<<" mass= "<<part.Mass()<<" p= "<<part.P()<<" E= "<<part.E()<<" trackID= "<<part.TrackId()<<" ND= "<<part.NumberDaughters()<<" Mother= "<<part.Mother()<<std::endl;
    
   
    if(part.PdgCode()==2212 && part.StatusCode()==1){
    //std::cout<<"and its mother = "<<part.Mother()<<std::endl;
    no_protons++;
    }
    
     if(part.PdgCode()==111 && part.StatusCode()==1){
    event_has_pi0=1;
    }
    
    if(part.PdgCode()==211 && part.StatusCode()==1){
    event_has_pi_plus=1;
    }
    
  
    
    
    }
 }
 
 std::cout<<"starting CCQE"<<std::endl;
 
  //..........................................................
  //          START	
  
  
  
  
  
 int no_mc_particles_in_event=0;
 double mu_theta_true, mu_phi_true;
 double pion_theta_true, pion_phi_true,momentum;
 double startpoint_x,startpoint_y, startpoint_z;
 double endpoint_x,endpoint_y, endpoint_z;
 std::string pri ("primary");
 double track_length=0;
 
 
 
 
 
 std::cout<<"geant_part.size()="<<geant_part.size()<<std::endl;
 
 for( unsigned int i = 0; i < geant_part.size(); ++i ){
   
  std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" trackId= "<<geant_part[i]->TrackId()<<" mother= "<<geant_part[i]->Mother()<<" process= "<<geant_part[i]->Process()<<std::endl;
     
     if(geant_part[i]->Process()==pri){std::cout<<" *****"<<std::endl;}

 momentum=sqrt(pow(geant_part[i]->Px(),2)+pow(geant_part[i]->Py(),2.)+pow(geant_part[i]->Pz(),2.));
 
    startpoint_x=geant_part[i]->Vx();
    startpoint_y=geant_part[i]->Vy();
    startpoint_z=geant_part[i]->Vz();
    
    endpoint_x=geant_part[i]->EndPoint()[0];
    endpoint_y=geant_part[i]->EndPoint()[1];
    endpoint_z=geant_part[i]->EndPoint()[2];
    track_length=sqrt(pow(endpoint_x-startpoint_x,2)+pow(endpoint_y-startpoint_y,2.)+pow(endpoint_z-startpoint_z,2.));
    std::cout<<"track_length= "<<track_length<<std::endl;

 //change the cut on energy?
 // count particles, exclude neutrons
     if(geant_part[i]->PdgCode()!= 2112 && geant_part[i]->Process()==pri && geant_part[i]->E()>0.01 && momentum> 0 && geant_part[i]->PdgCode()<10000){
  
  proton_track_length->Fill(track_length);
  
  if(track_length>1){no_mc_particles_in_event++;}
  
   std::cout<<"counting as a valid particle pdg= "<<geant_part[i]->PdgCode()<< "no_mc_particles_in_event= "<<no_mc_particles_in_event<<std::endl;
   
    }
   
   momentum=sqrt(pow(geant_part[i]->Px(),2)+pow(geant_part[i]->Py(),2.)+pow(geant_part[i]->Pz(),2.));
    
    //WORK ON PROTONS:
//..........................................................................
    if((geant_part[i]->PdgCode()==2212) && (geant_part[i]->Process()==pri) && momentum >0 && track_length>1){
    std::cout<<"we have a proton with E="<<geant_part[i]->E()<<" and mom= "<<momentum<<std::endl;
    have_p++;
    E_p=geant_part[i]->E();
   
    
    
    //check if protons are contained in the detector (check their exit point):
    std::cout<<"vertex [x,y,z]= [ "<<geant_part[i]->Vx()<<", "<<geant_part[i]->Vy()<<", "<<geant_part[i]->Vz()<<" ]" <<" Exit point of proton [x,y,z]= [ "<<geant_part[i]->EndPoint()[0]<<", "<<geant_part[i]->EndPoint()[1]<<", "<<geant_part[i]->EndPoint()[2]<<" ]"<<std::endl;
    
    if(endpoint_x >3 && endpoint_x <44 && endpoint_y >-17 && endpoint_y <17 && endpoint_z >3 && endpoint_z <87){
     std::cout<<"proton fully contained"<<std::endl;
     P_containment->Fill(1);
    }
    else{
    std::cout<<"proton NOT contained!!!"<<std::endl;
     P_containment->Fill(0);
    }
    
    
    }
  //..........................................................................
  
    if((geant_part[i]->PdgCode()==13 || geant_part[i]->PdgCode()==-13) && (geant_part[i]->Process()==pri)){
   
   
   if(geant_part[i]->PdgCode()==13){
   mu_charge->Fill(-1);
   }
    if(geant_part[i]->PdgCode()==-13){
   mu_charge->Fill(1);
   }
   
   
    std::cout<<"we have a mu with E="<<geant_part[i]->E()<<std::endl;
    
    mu_theta_true=(TMath::ACos(geant_part[i]->Pz()/sqrt(pow(geant_part[i]->Px(),2)+pow(geant_part[i]->Py(),2)+pow(geant_part[i]->Pz(),2))))*(180/TMath::Pi());
    std::cout<<"mu_theta_true= "<<mu_theta_true<<std::endl;
    
    mu_phi_true=(TMath::Pi()+TMath::ATan2(-geant_part[i]->Py(),-geant_part[i]->Px()))*(180/TMath::Pi());
 std::cout<<"mu_phi_true= "<<mu_phi_true<<std::endl;
 
 Mu_theta->Fill(mu_theta_true);
 Mu_phi->Fill(mu_phi_true);
 // if(mu_phi_true<=180){
//  Mu_phi_oneside->Fill(mu_phi_true);}
//  if(mu_phi_true>180){
//  Mu_phi_oneside->Fill(360-mu_phi_true);}
 
 
    }
    
   
    
 
 
 }
 
No_protons_in_event->Fill(have_p);

 No_particles_in_event->Fill(no_mc_particles_in_event);
    


have_p=0;
no_mc_particles_in_event=0;
geant_part.clear();
mclist.clear();








 
}

