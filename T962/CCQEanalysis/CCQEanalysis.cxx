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

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "T962/CCQEanalysis/CCQEanalysis.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"


 
//-------------------------------------------------
t962::CCQEanalysis::CCQEanalysis(fhicl::ParameterSet const& pset) : 
 
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fClusterFinderModuleLabel (pset.get< std::string >("ClusterFinderModuleLabel")),
  fDBClusterFinderModuleLabel (pset.get< std::string >("DBClusterFinderModuleLabel")),
  fHoughModuleLabel (pset.get< std::string >("HoughModuleLabel")),
  fLineMModuleLabel       (pset.get< std::string >("LineMModuleLabel")      ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  //fMINOSModuleLabel         (pset.get< std::string >("MINOSModuleLabel")        ),
  //fTrackMatchModuleLabel    (pset.get< std::string >("TrackMatchModuleLabel")  
  frun(0),
  fevent(0),
  fbeam(0),
  fno_primaries(100),
  no_hits(5000)
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
 delete hit_plane;
 delete hit_wire;
 delete hit_channel;
 delete hit_peakT;
 delete hit_charge;
 delete twodvtx_w_reco;
delete twodvtx_t_reco;
 delete twodvtx_w_truth;
delete twodvtx_t_truth;
delete fNumberDaughters;

 
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
  
  hit_plane= new int[no_hits];
   hit_wire= new int[no_hits];
    hit_channel= new int[no_hits];
   hit_peakT= new double[no_hits];
   hit_charge= new double[no_hits];
  twodvtx_w_reco= new double[2];
  twodvtx_t_reco= new double[2];
  twodvtx_w_truth= new double[2];
  twodvtx_t_truth= new double[2];
  fNumberDaughters= new int[fno_primaries];
  
  
  
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
  
  fTree->Branch("no_hits",&no_hits,"no_hits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[no_hits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[no_hits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[no_hits]/I");
   fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/D");
   fTree->Branch("hit_charge",hit_charge,"hit_charge[no_hits]/D");
//   
  fTree->Branch("twodvtx_w_reco", twodvtx_w_reco, "twodvtx_w_reco[2]/D");
  fTree->Branch("twodvtx_t_reco", twodvtx_t_reco, "twodvtx_t_reco[2]/D");
  fTree->Branch("twodvtx_w_truth", twodvtx_w_truth, "twodvtx_w_truth[2]/D");
  fTree->Branch("twodvtx_t_truth", twodvtx_t_truth, "twodvtx_t_truth[2]/D");
  
  fTree->Branch("NumberDaughters",fNumberDaughters,"NumberDaughters[no_primaries]/I");
  
//                  RECO :
  fTree->Branch("vtxx_reco",&vtxx_reco,"vtxx_reco/D");
  fTree->Branch("vtxy_reco",&vtxy_reco,"vtxy_reco/D");
  fTree->Branch("vtxz_reco",&vtxz_reco,"vtxz_reco/D");
  fTree->Branch("ndbclusu_reco",&ndbclusu_reco,"ndbclusu_reco/I");
  fTree->Branch("ndbclusv_reco",&ndbclusv_reco,"ndbclusv_reco/I");
  fTree->Branch("ndbclusw_reco",&ndbclusw_reco,"ndbclusw_reco/I");
  fTree->Branch("nhoughu_reco",&nhoughu_reco,"nhoughu_reco/I");
  fTree->Branch("nhoughv_reco",&nhoughu_reco,"nhoughv_reco/I");
  fTree->Branch("nhoughw_reco",&nhoughu_reco,"nhoughw_reco/I");
  fTree->Branch("nlineu_reco",&nlineu_reco,"nlineu_reco/I");
  fTree->Branch("nlinev_reco",&nlinev_reco,"nlinev_reco/I");
  fTree->Branch("nlinew_reco",&nlinew_reco,"nlinew_reco/I");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("nclusu_reco",&nclusu_reco,"nclusu_reco/I");
  fTree->Branch("nclusv_reco",&nclusv_reco,"nclusv_reco/I");
  fTree->Branch("nclusw_reco",&nclusw_reco,"nclusw_reco/I");
  // fTree->Branch("nvertextracks_reco",&nvertextracks_reco,"nvertextracks_reco/I");
//   fTree->Branch("nvertexclustersu_reco",&nvertexclustersu_reco,"nvertexclustersu_reco/I");
//   fTree->Branch("nvertexclustersv_reco",&nvertexclustersv_reco,"nvertexclustersv_reco/I");
//   fTree->Branch("nvertexclustersw_reco",&nvertexclustersw_reco,"nvertexclustersw_reco/I");
// 
//  fTree->Branch("trackstart_dcosx_reco",&trackstart_dcosx_reco, "trackstart_dcosx_reco/D");
//   fTree->Branch("trackstart_dcosy_reco",&trackstart_dcosy_reco, "trackstart_dcosy_reco/D");
//   fTree->Branch("trackstart_dcosz_reco",&trackstart_dcosz_reco, "trackstart_dcosz_reco/D");
//   fTree->Branch("trackexit_dcosx_reco",&trackexit_dcosx_reco, "trackexit_dcosx_reco/D");
//   fTree->Branch("trackexit_dcosy_reco",&trackexit_dcosy_reco, "trackexit_dcosy_reco/D");
//   fTree->Branch("trackexit_dcosz_reco",&trackexit_dcosz_reco, "trackexit_dcosz_reco/D");
//   fTree->Branch("trackstart_x_reco",&trackstart_x_reco, "trackstart_x_reco/D");
//   fTree->Branch("trackstart_y_reco",&trackstart_y_reco, "trackstart_y_reco/D");
//   fTree->Branch("trackstart_z_reco",&trackstart_z_reco, "trackstart_z_reco/D");
//   fTree->Branch("trackexit_x_reco",&trackexit_x_reco, "trackexit_x_reco/D");
//   fTree->Branch("trackexit_y_reco",&trackexit_y_reco, "trackexit_y_reco/D");
//   fTree->Branch("trackexit_z_reco",&trackexit_z_reco, "trackexit_z_reco/D");    
//   fTree->Branch("nmatched_reco",&nmatched_reco,"nmatched_reco/I");  
//   fTree->Branch("trk_mom_minos",&trk_mom_minos,"trk_mom_minos/D");
//   fTree->Branch("trk_charge_minos",&trk_charge_minos,"trk_charge_minos/D");
//   fTree->Branch("trk_dcosx_minos",&trk_dcosx_minos,"trk_dcosx_minos/D");
//   fTree->Branch("trk_dcosy_minos",&trk_dcosy_minos,"trk_dcosy_minos/D");
//   fTree->Branch("trk_dcosz_minos",&trk_dcosz_minos,"trk_dcosz_minos/D");
//   fTree->Branch("trk_vtxx_minos",&trk_vtxx_minos,"trk_vtxx_minos/D");
//   fTree->Branch("trk_vtxy_minos",&trk_vtxy_minos,"trk_vtxy_minos/D");
//   fTree->Branch("trk_vtxz_minos",&trk_vtxz_minos,"trk_vtxz_minos/D");
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
   fbeam=-1; //1 for nu mode, -1 for anti-nu mode
   
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

  sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArG4ModuleLabel);
  
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  std::vector< art::Ptr<recob::Hit> > hits;
  art::fill_ptr_vector(hits, hitListHandle);
  
  art::Handle< std::vector<recob::Cluster> > dbclusterListHandle;
  evt.getByLabel(fDBClusterFinderModuleLabel,dbclusterListHandle);
 
   art::Handle< std::vector<sim::Particle> > geant_list;
    evt.getByLabel (fLArG4ModuleLabel,geant_list);
    
  //.................
  
   art::Handle< std::vector<recob::Cluster> > linemListHandle;
  evt.getByLabel(fLineMModuleLabel,linemListHandle);
  
   art::Handle< std::vector<recob::Cluster> > houghListHandle;
  evt.getByLabel(fHoughModuleLabel,houghListHandle);
  
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  //art::Handle< std::vector<t962::MINOS> > minosListHandle;
  //evt.getByLabel(fMINOSModuleLabel,minosListHandle);
  //art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
  //evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);
  

  
  //---------------------------------------------------
art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
  
  art::PtrVector<recob::Cluster> clusters;
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }
    //---------------------------------------------------
 std::cout<<"-----------------------------------"<<std::endl;
  
   nclusu_reco=0;
   nclusv_reco=0;
   nclusw_reco=0;
  
  
  std::cout<<"*** DBSCAN list size= "<<clusters.size()<<std::endl;
  for(unsigned int i=0; i<clusters.size();++i){
   switch(clusters[i]->View()){
    case geo::kU :
      nclusu_reco ++;
     // Cls[0].push_back(i);
      std::cout<<"here1"<<std::endl;
      break;
    case geo::kV :
      nclusv_reco ++;
      //Cls[1].push_back(i);
      std::cout<<"here2"<<std::endl;
      break;
    case geo::kW :
      nclusw_reco ++;
      //Cls[2].push_back(i);
      std::cout<<"here3"<<std::endl;
      break;
    default :
      break;
    }
  }
std::cout<<"No dbscan in u= "<<nclusu_reco<<std::endl;
std::cout<<"No dbscan in v= "<<nclusv_reco<<std::endl;
std::cout<<"No dbscan in w= "<<nclusw_reco<<std::endl;
  
  
   std::cout<<"-----------------------------------"<<std::endl;
   
   
   
  
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
    
    
  art::PtrVector<recob::Cluster> dbclusterlist;
  dbclusterlist.clear();
  if(evt.getByLabel(fDBClusterFinderModuleLabel,dbclusterListHandle))
  for (unsigned int ii = 0; ii <  dbclusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> dbclusterHolder(linemListHandle,ii);
      dbclusterlist.push_back(dbclusterHolder);
    }
    
    art::PtrVector<recob::Cluster> houghlist;
  if(evt.getByLabel(fHoughModuleLabel,houghListHandle))
  for (unsigned int ii = 0; ii <  houghListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> houghHolder(houghListHandle,ii);
      houghlist.push_back(houghHolder);
    }
    
    
    art::PtrVector<recob::Cluster> linemlist;
  if(evt.getByLabel(fLineMModuleLabel,linemListHandle))
  for (unsigned int ii = 0; ii <  linemListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> linemHolder(linemListHandle,ii);
      linemlist.push_back(linemHolder);
    }

  art::PtrVector<recob::Track> tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
  for (unsigned int i = 0; i < trackListHandle->size(); ++i){
    art::Ptr<recob::Track> trackHolder(trackListHandle,i);
    tracklist.push_back(trackHolder);
  }

  art::PtrVector<recob::EndPoint2D> endpointlist;
  if(evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle))
    for (unsigned int i = 0; i < endpointListHandle->size(); ++i){
      art::Ptr<recob::EndPoint2D> endpointHolder(endpointListHandle,i);
      endpointlist.push_back(endpointHolder);
    }

  art::PtrVector<recob::Vertex> vertexlist;
  if(evt.getByLabel(fVertexModuleLabel,vertexListHandle))
  for (unsigned int i = 0; i < vertexListHandle->size(); ++i){
    art::Ptr<recob::Vertex> vertexHolder(vertexListHandle,i);
    vertexlist.push_back(vertexHolder);
  }

 //  art::PtrVector<t962::MINOS> minoslist;
//   if(evt.getByLabel(fMINOSModuleLabel,minosListHandle))
//   for (unsigned int i = 0; i < minosListHandle->size(); i++){
//     art::Ptr<t962::MINOS> minosHolder(minosListHandle,i);
//     minoslist.push_back(minosHolder);
//   }
// 
//   art::PtrVector<t962::MINOSTrackMatch> trackmatchlist;
//   if(evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle))
//   for (unsigned int i = 0; i < trackmatchListHandle->size(); i++){
//     art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,i);
//     trackmatchlist.push_back(trackmatchHolder);
//   }
  
  
     // 2d vertex information
  bool found2dvtx = false;
  
  for (unsigned int j = 0; j<endpointlist.size();j++){
 std::cout<<"j="<<j<<" W_VERTEX_RECO= "<<endpointlist[j]->WireNum()<<" T_VERTEX_RECO= "<<endpointlist[j]->DriftTime()<<std::endl;
          twodvtx_w_reco[j]=endpointlist[j]->WireNum();
          twodvtx_t_reco[j]=endpointlist[j]->DriftTime();
          found2dvtx=true;
    }
  
  
   //................................................................
   art::ServiceHandle<geo::Geometry> geom;  
   
double vertex [3] = { 0, 0, 0 };
int have_p=0;
int have_pion=0;
double E_p, E_pion;
int event_has_pi0=0;
int event_has_pi_plus=0;

double MC_Total_Eng=0;

  
    
    
    
    
  //..........................................................
   
    //First determine what kind of event we are dealing with, and select only the ones you want:
 //..........................................................
 int no_protons=0;
 int mode=-999;
 int ccnc=-999;
 
 
 for( unsigned int i = 0; i < 1; ++i ){
 //for( unsigned int i = 0; i < mclist.size(); ++i ){
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
 
 
 
 
 //.........................................................
 
   //get true 2d vertex:
     art::ServiceHandle<util::LArProperties> larp;
    
    for( unsigned int i = 0; i < 1; ++i ){
    //for( unsigned int i = 0; i < mclist.size(); ++i ){

    art::Ptr<simb::MCTruth> mc(mclist[i]);

    simb::MCParticle neut(mc->GetParticle(i));

    
    fMCvertex[0] =neut.Vx();
    fMCvertex[1] =neut.Vy();
    fMCvertex[2] =neut.Vz();
   
    double presamplings=60.0;
    double drifttick=(fMCvertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198)+presamplings;
    
    twodvtx_t_truth[0]=drifttick;
    twodvtx_t_truth[1]=drifttick;
  }
  // now wire vertex:
   unsigned int channel2,plane2,wire2,tpc2; 
  for(size_t tpc = 0; tpc < geom->NTPC(); ++tpc){
   
  for(unsigned int plane=0;plane<geom->Nplanes(tpc);plane++){
  if(plane==0){
	fMCvertex[0]=.3;//force time coordinate to be closer to induction plane 
	}
      else{
	fMCvertex[0]=-.3;//force time coordinate to be closer to collection plane
     }
     
  
  try{
  channel2 = geom->NearestChannel(fMCvertex,plane);
  }
  catch(cet::exception &e){
  mf::LogWarning("KingaClusterexc")<<e;
  std::cout<<std::endl;
   std::cout<<std::endl;
    std::cout<<std::endl;
     std::cout<<"exception *****************************"<<std::endl;
     std::cout<<std::endl;
   std::cout<<std::endl;
    std::cout<<std::endl;
    
  std::cout<<"fMCvertex[2]= "<<fMCvertex[2]<<" plane="<<plane<<std::endl;
  if(plane==0 && fMCvertex[2]<5) channel2=0;
  else if(plane==0 && fMCvertex[2]>geom->DetLength()-5) channel2=(geom->Nchannels())/2 -1;
  else if(plane==1 && fMCvertex[2]>geom->DetLength()-5) channel2=geom->Nchannels()-1;
  else if(plane==1 && fMCvertex[2]<5) channel2=(geom->Nchannels())/2 -1;

  
  }
      geom->ChannelToWire(channel2,tpc2,plane2,wire2);   
   
   
   twodvtx_w_truth[plane]=wire2;
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
   
    fprimaries_pdg[i]=geant_part[i]->PdgCode();
    
    fEng[i]=geant_part[i]->E();
    fPx[i]=geant_part[i]->Px();
   
    fPy[i]=geant_part[i]->Py();
    fPz[i]=geant_part[i]->Pz();
    
   fStartPointx[i]=geant_part[i]->Vx();
   fStartPointy[i]=geant_part[i]->Vy();
   fStartPointz[i]=geant_part[i]->Vz();
   fEndPointx[i]=geant_part[i]->EndPoint()[0];
   fEndPointy[i]=geant_part[i]->EndPoint()[1];
   fEndPointz[i]=geant_part[i]->EndPoint()[2];
   
   fNumberDaughters[i]=geant_part[i]->NumberDaughters();
   
   std::cout<<"geant_part[i]->EndPoint()[0]= "<<geant_part[i]->EndPoint()[0]<<std::endl;
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
 
 
 
 //.....................................................................
  
 //                RECONSTRUCTION INFO:
 //.....................................................................

 
 //vertex information
  if(vertexlist.size())
  {
    double vtxxyz[3];
    vertexlist[0]->XYZ(vtxxyz);
    vtxx_reco = vtxxyz[0];
    vtxy_reco = vtxxyz[1];
    vtxz_reco = vtxxyz[2];
  }
  // DBSCANcluster information
  ndbclusu_reco = 0;
  ndbclusv_reco = 0;
  ndbclusw_reco = 0;

  int nplanes = geom->Nplanes();
  std::vector<int> Cls[nplanes];
std::cout<<"*** DBSCAN list size= "<<dbclusterlist.size()<<std::endl;
  //for(unsigned int i=0; i<dbclusterlist.size();++i)// {
//   // std::cout<<"i= "<<i<<std::endl;
// //   if(clusterlist[i]->View()==geo::kU){
// //   std::cout<<" in U?"<<std::endl;
// //   nclusu_reco++;
// //   }
// //   else if( clusterlist[i]->View()==geo::kV){
// //   std::cout<<" in V?"<<std::endl;
// //   nclusv_reco++;
// //   }
// //   else{
// //   
// //   std::cout<<"no idea <<<<<<<<<<<<"<<std::endl;
// //   }
// 
// 
//    switch(dbclusterlist[i]->View()){
//     case geo::kU :
//       ndbclusu_reco ++;
//       Cls[0].push_back(i);
//       std::cout<<"here1"<<std::endl;
//       break;
//     case geo::kV :
//       ndbclusv_reco ++;
//       Cls[1].push_back(i);
//       std::cout<<"here2"<<std::endl;
//       break;
//     case geo::kW :
//       ndbclusw_reco ++;
//       Cls[2].push_back(i);
//       std::cout<<"here3"<<std::endl;
//       break;
//     default :
//       break;
//     }
//   }
// std::cout<<"No dbscan in u= "<<ndbclusu_reco<<std::endl;
// std::cout<<"No dbscan in v= "<<ndbclusv_reco<<std::endl;
// std::cout<<"No dbscan in w= "<<ndbclusw_reco<<std::endl;

// HOUGH cluster information
  nhoughu_reco = 0;
  nhoughv_reco = 0;
  nhoughw_reco = 0;

 
  std::vector<int> Hough[nplanes];

  for(unsigned int i=0; i<houghlist.size();++i){
    
    switch(houghlist[i]->View()){
    case geo::kU :
      nhoughu_reco ++;
      Hough[0].push_back(i);
      break;
    case geo::kV :
      nhoughv_reco ++;
      Hough[1].push_back(i);
      break;
    case geo::kW :
      nhoughw_reco ++;
      Hough[2].push_back(i);
      break;
    default :
      break;
    }
  }

std::cout<<"No hough in u= "<<nhoughu_reco<<std::endl;
std::cout<<"No hough in v= "<<nhoughv_reco<<std::endl;
std::cout<<"No hough in w= "<<nhoughw_reco<<std::endl;

 
  // LINE MERGER cluster information
  nlineu_reco = 0;
  nlinev_reco = 0;
  nlinew_reco = 0;

  
  std::vector<int> Lines[nplanes];

  for(unsigned int i=0; i<linemlist.size();++i){
    
    switch(linemlist[i]->View()){
    case geo::kU :
      nlineu_reco ++;
      Lines[0].push_back(i);
      break;
    case geo::kV :
      nlinev_reco ++;
      Lines[1].push_back(i);
      break;
    case geo::kW :
      nlinew_reco ++;
      Lines[2].push_back(i);
      break;
    default :
      break;
    }
  }
 
 std::cout<<"No line mergers in u= "<<nlineu_reco<<std::endl;
 std::cout<<"No line mergers in v= "<<nlinev_reco<<std::endl;
 std::cout<<"No line mergers in w= "<<nlinew_reco<<std::endl;
 
 
 
 
 //track information
     ntracks_reco=tracklist.size();
     std::cout<<"No of 3d tracks= "<<ntracks_reco<<std::endl;
 
 
//////////////////////////////////////////////////////////////////////////
/////////        HIT info:       ////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
  std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
   unsigned int p(0),w(0),t(0), channel(0);
   int hit_no;
   
   
   std::cout<<"hits.size()= "<<hits.size()<<std::endl;
   
   no_hits=hits.size();
   
  while( itr != hits.end() ){
 // std::cout<<"working on hit# "<<itr-hits.begin()<<std::endl;
  hit_no=(int)(itr-hits.begin());
  //std::cout<<"hit_no= "<<hit_no<<std::endl;
  //std::cout<<"channel= "<<(*itr)->Channel()<<" ";
  channel=(*itr)->Channel();
  
   geom->ChannelToWire(channel,t,p,w);
   hit_channel[hit_no]= channel;
  hit_plane[hit_no]=p;
   hit_wire[hit_no]=w;
   hit_peakT[hit_no]=(*itr)->PeakTime();
   hit_charge[hit_no]=(*itr)->Charge();
   
   //std::cout<<"p= "<<p<<" ,w= "<<w<<" ,time= "<<(*itr)->PeakTime()<<" Charge= "<<(*itr)->Charge()<<std::endl;
   
  //  if((*itr)->View()==geo::kU){
//     std::cout<<"plane 0"<<std::endl;}
//     else if((*itr)->View()==geo::kV){
//     std::cout<<"plane 1"<<std::endl;}
  itr++;
  }
 
 
 ////////////////////////////////////////////////////////////////////////
 fTree->Fill();
//dbclusterlist.clear();
  
 
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
