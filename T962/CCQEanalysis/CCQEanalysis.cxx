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

#include "TTree.h"

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
  fClusterFinderModuleLabel (pset.get< std::string >("ClusterFinderModuleLabel")),
  frun(0),
  fevent(0),
  fbeam(0),
  fno_primaries(100)
{



}

//-------------------------------------------------
t962::CCQEanalysis::~CCQEanalysis()
{
delete fprimaries_pdg;
delete fEng;
delete fPx;
 delete fPy;
 delete fPz;
 delete fStartPointx;
 delete fStartPointy;
 delete fStartPointz;
 delete fEndPointx;
 delete fEndPointy;
 delete fEndPointz;
 
}

void t962::CCQEanalysis::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  
  //..........................................................
  
  fTree = tfs->make<TTree>("anatree","ccqe analysis");
  fprimaries_pdg= new int[fno_primaries];
  fEng= new double[fno_primaries];
  fPx= new double[fno_primaries];
  fPy= new double[fno_primaries];
  fPz= new double[fno_primaries];
  fStartPointx= new double[fno_primaries];
  fStartPointy= new double[fno_primaries];
  fStartPointz= new double[fno_primaries];
  fEndPointx= new double[fno_primaries];
  fEndPointy= new double[fno_primaries];
  fEndPointz= new double[fno_primaries];
  
   fTree->Branch("beam",&fbeam,"beam/I");
   fTree->Branch("run",&frun,"run/I");
   fTree->Branch("event",&fevent,"event/I");
  //......................................................
  //fTree->Branch("pot",&pot,"pot/D");
  //fTree->Branch("isdata",&isdata,"isdata/I");
  // fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
//   fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
//   fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
//   fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
//   fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
//   fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");
//   
  //......................................................

   fTree->Branch("no_primaries",&fno_primaries,"no_primaries/I");
   
   
   //brian:
  // fTree->Branch("primaries_pdg",&fprimaries_pdg);
//   fTree->Branch("Eng",&fEng);
  //...................................................
  
  
  fTree->Branch("primaries_pdg",fprimaries_pdg,"primaries_pdg[no_primaries]/I");
  fTree->Branch("Eng",fEng,"Eng[no_primaries]/D");
  
   fTree->Branch("Px",fPx,"Px[no_primaries]/D");
 fTree->Branch("Py",fPy,"Py[no_primaries]/D");
  fTree->Branch("Pz",fPz,"Pz[no_primaries]/D");
  fTree->Branch("StartPointx",fStartPointx,"StartPointx[no_primaries]/D");
  fTree->Branch("StartPointy",fStartPointy,"StartPointy[no_primaries]/D");
  fTree->Branch("StartPointz",fStartPointz,"StartPointz[no_primaries]/D");
  fTree->Branch("EndPointx",fEndPointx,"EndPointx[no_primaries]/D");
  fTree->Branch("EndPointy",fEndPointy,"EndPointy[no_primaries]/D");
  fTree->Branch("EndPointz",fEndPointz,"EndPointz[no_primaries]/D");
  
//   
  
  
  // fTree->Branch("theta_truth",&Mu_theta_truth,"Mu_theta_truth/D");
//   fTree->Branch("phi_truth",&Mu_phi_truth,"Mu_phi_truth/D");
//   fTree->Branch("charge_truth",&Mu_charge_truth,"Mu_charge_truth/I");
//   
  
  
  
  
  
  
  
    //..........................................................

Mu_theta=tfs->make<TH1F>("Mu_theta","Muon theta angle", 360,0 ,360);
Mu_phi=tfs->make<TH1F>("Mu_phi","Muon phi angle", 360,0 ,360);
mu_charge=tfs->make<TH1F>("mu_charge","Muon charge", 4,-2 ,2);
No_protons_in_event=tfs->make<TH1F>("No_protons_in_event","No of protons in each CCQE event", 10,0 ,10);
No_particles_in_event=tfs->make<TH1F>("No_particles_in_event","No of particles in each CCQE event (excluding neutrons), track>1 cm", 20,0 ,20);
P_containment=tfs->make<TH1F>("P_containment","Proton contained=1, else 0", 3,-1 ,2);
P_containment_no_cut=tfs->make<TH1F>("P_containment_no_cut","Proton contained=1, else 0, no active volume cut", 3,-1 ,2);
proton_track_length=tfs->make<TH1F>("proton_track_length","Length of proton track", 300,0 ,150);

 track_l_vs_containment = tfs->make<TH1F>("track_l_vs_containment", ";Containment; Proton Track Length", 4,  -1, 3);

 p_mom_vs_containment = tfs->make<TH1F>("p_mom_vs_containment", ";Containment; Proton Momentum", 4,  -1, 3);
  
  track_l_vs_mom = tfs->make<TH1F>("track_l_vs_mom", ";Proton Track Length; Momentum",5000, 0, 5);

Vertex_x=tfs->make<TH1F>("Vertex_x","vertex X coordinate value", 70,-10 ,60);
Vertex_y=tfs->make<TH1F>("Vertex_y","vertex Y coordinate value", 60,-30 ,30);
Vertex_z=tfs->make<TH1F>("Vertex_z","vertex Z coordinate value", 130,-20 ,110);

  
}

void t962::CCQEanalysis::analyze(const art::Event& evt)
{

   frun = evt.run();
   fevent = evt.id().event();
   fbeam=1; //1 for nu mode, -1 for anti-nu mode
   
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
 int mode=-999;
 int ccnc=-999;
 
 
 for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
    simb::MCParticle neut(mc->GetParticle(i));
    
    ccnc=mclist[0]->GetNeutrino().CCNC();
    mode=mclist[0]->GetNeutrino().Mode();
    std::cout<<"ccnc= "<<ccnc<<std::endl;
    std::cout<<"Mode= "<<mode<<std::endl;
    
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
 
 
 
  //..........................................................
  //          START	
  
 int primary=0;
 int no_mc_particles_in_event=0;
 double mu_theta_true, mu_phi_true;
 double pion_theta_true, pion_phi_true,momentum;
 double startpoint_x,startpoint_y, startpoint_z;
 double endpoint_x,endpoint_y, endpoint_z;
 std::string pri ("primary");
 double track_length=0;
 
  //determine the number of primary particles from geant:
  
  for( unsigned int i = 0; i < geant_part.size(); ++i ){
   
    if(geant_part[i]->Process()==pri){
    primary++;
    }
   
   }
  
  fno_primaries=primary;
  
  

 
 
 // fprimaries_pdg.clear();
//  fEng.clear();
 
 std::cout<<"geant_part.size()="<<geant_part.size()<<std::endl;
 
 for( unsigned int i = 0; i < geant_part.size(); ++i ){
   
   if(geant_part[i]->Process()==pri){
   
   
    // fprimaries_pdg.push_back(geant_part[i]->PdgCode());
//     std::cout<<"geant_part[i]->E()= "<<geant_part[i]->E()<<std::endl;
//     fEng.push_back(geant_part[i]->E());
//     
    std::cout<<"here1"<<std::endl;
    fprimaries_pdg[i]=geant_part[i]->PdgCode();
    std::cout<<"here2"<<std::endl;
    fEng[i]=geant_part[i]->E();
    fPx[i]=geant_part[i]->Px();
    std::cout<<"here3"<<std::endl;
    fPy[i]=geant_part[i]->Py();
    fPz[i]=geant_part[i]->Pz();
    
   fStartPointx[i]=geant_part[i]->Vx();
   fStartPointy[i]=geant_part[i]->Vy();
   fStartPointz[i]=geant_part[i]->Vz();
   fEndPointx[i]=geant_part[i]->EndPoint()[0];
   fEndPointy[i]=geant_part[i]->EndPoint()[1];
   fEndPointz[i]=geant_part[i]->EndPoint()[2];
   
  std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" trackId= "<<geant_part[i]->TrackId()<<" mother= "<<geant_part[i]->Mother()<<" process= "<<geant_part[i]->Process()<<std::endl;
     
     

 momentum=sqrt(pow(geant_part[i]->Px(),2)+pow(geant_part[i]->Py(),2.)+pow(geant_part[i]->Pz(),2.));
 std::cout<<"*** mom= "<<momentum<<std::endl;
 
    startpoint_x=geant_part[i]->Vx();
    startpoint_y=geant_part[i]->Vy();
    startpoint_z=geant_part[i]->Vz();
    
    endpoint_x=geant_part[i]->EndPoint()[0];
    endpoint_y=geant_part[i]->EndPoint()[1];
    endpoint_z=geant_part[i]->EndPoint()[2];
    track_length=sqrt(pow(endpoint_x-startpoint_x,2)+pow(endpoint_y-startpoint_y,2.)+pow(endpoint_z-startpoint_z,2.));
    std::cout<<"track_length= "<<track_length<<std::endl;

 
 // count particles, exclude neutrons
     if(geant_part[i]->PdgCode()!= 2112 && geant_part[i]->Process()==pri && momentum> 0 && geant_part[i]->PdgCode()<10000){
  
  if(geant_part[i]->PdgCode()==2212){
  proton_track_length->Fill(track_length);
  track_l_vs_mom->Fill(track_length,momentum);
  }
  std::cout<<"before, no_mc_particles_in_event= "<<no_mc_particles_in_event<<std::endl;
  
  if(track_length>1.5){no_mc_particles_in_event++;}
  
   std::cout<<"counting as a valid particle pdg= "<<geant_part[i]->PdgCode()<< " no_mc_particles_in_event= "<<no_mc_particles_in_event<<std::endl;
   
    }
    
   //why is this written again?:
   //momentum=sqrt(pow(geant_part[i]->Px(),2)+pow(geant_part[i]->Py(),2.)+pow(geant_part[i]->Pz(),2.));
    
   
    //WORK ON PROTONS:
//..........................................................................
    if((geant_part[i]->PdgCode()==2212) && (geant_part[i]->Process()==pri) && momentum >0 && track_length>1.5){
    std::cout<<"we have a proton with E="<<geant_part[i]->E()<<" and mom= "<<momentum<<std::endl;
    have_p++;
    E_p=geant_part[i]->E();
   
    
    
    //check if protons are contained in the detector (check their exit point):
    std::cout<<"vertex [x,y,z]= [ "<<geant_part[i]->Vx()<<", "<<geant_part[i]->Vy()<<", "<<geant_part[i]->Vz()<<" ]" <<" Exit point of proton [x,y,z]= [ "<<geant_part[i]->EndPoint()[0]<<", "<<geant_part[i]->EndPoint()[1]<<", "<<geant_part[i]->EndPoint()[2]<<" ]"<<std::endl;
    
    if(endpoint_x >3 && endpoint_x <44 && endpoint_y >-17 && endpoint_y <17 && endpoint_z >3 && endpoint_z <87){
     std::cout<<"proton fully contained"<<std::endl;
     P_containment->Fill(1);
     track_l_vs_containment->Fill(1.0,track_length);
     p_mom_vs_containment->Fill(1.0,momentum);
     
    }
    else{
    std::cout<<"proton NOT contained!!!"<<std::endl;
     P_containment->Fill(0);
     track_l_vs_containment->Fill(0.,track_length);
     p_mom_vs_containment->Fill(0.,momentum);
    }
    
     if(endpoint_x >0 && endpoint_x <47 && endpoint_y >-20 && endpoint_y <20 && endpoint_z >0 && endpoint_z <90){
     std::cout<<"proton fully contained"<<std::endl;
     P_containment_no_cut->Fill(1);
     
    }
    
    else{
    std::cout<<"proton NOT contained!!!"<<std::endl;
     P_containment_no_cut->Fill(0);
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
    
   
     
  }//if primary process
 
 } //loop thru geant particles
 
 fTree->Fill();

  
 
No_protons_in_event->Fill(have_p);

 No_particles_in_event->Fill(no_mc_particles_in_event);
    
 


// ccnc_truth = mclist[0]->GetNeutrino().CCNC();
// mode_truth = mclist[0]->GetNeutrino().Mode();
// 
//  if (mclist[0]->GetNeutrino().Lepton().P()){
//  lep_mom_truth = mclist[0]->GetNeutrino().Lepton().P();
//       lep_dcosx_truth = mclist[0]->GetNeutrino().Lepton().Px()/mclist[0]->GetNeutrino().Lepton().P();
//       lep_dcosy_truth = mclist[0]->GetNeutrino().Lepton().Py()/mclist[0]->GetNeutrino().Lepton().P();
//       lep_dcosz_truth = mclist[0]->GetNeutrino().Lepton().Pz()/mclist[0]->GetNeutrino().Lepton().P();
//  
//  
//  }
 

//fTree->Fill();

have_p=0;
no_mc_particles_in_event=0;
//geant_part.clear();
//mclist.clear();



 
}

// void t962::CCQEanalysis::ResetVars(){
// 
//   run = -99999;
//   event = -99999;
//   ccnc_truth = -99999;
//   mode_truth = -99999;
//   lep_mom_truth = -99999;
//   lep_dcosx_truth = -99999;
//   lep_dcosy_truth = -99999;
//   lep_dcosz_truth = -99999;
//   
//   
//   }