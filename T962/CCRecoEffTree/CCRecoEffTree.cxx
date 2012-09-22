////////////////////////////////////////////////////////////////////////
//
// Create a TTree for CC nu-mu Event Reconstruction Efficiency Study
//
// \author saima@ksu.edu
// 
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>

#include "TTree.h"

#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/AssociationUtil.h" 
#include "art/Framework/Core/FindMany.h"


#include "T962/CCRecoEffTree/CCRecoEffTree.h"
#include "T962/T962_Objects/MINOS.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "SummaryData/POTSummary.h"
#include "MCCheater/BackTracker.h"

 
//-------------------------------------------------
t962::CCRecoEffTree::CCRecoEffTree(fhicl::ParameterSet const& pset) : 
  
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fMINOSModuleLabel         (pset.get< std::string >("MINOSModuleLabel")        ),
  fTrackMatchModuleLabel    (pset.get< std::string >("TrackMatchModuleLabel")   ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fboundaryWindow           (pset.get< double >("boundaryWindow")               )
  
{
}

//-------------------------------------------------
t962::CCRecoEffTree::~CCRecoEffTree()
{
  delete trk_charge_minos_all; 
}

void t962::CCRecoEffTree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("EffTree","EffTree");

  trk_charge_minos_all = new double[nminos_tracks];

  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("vtxx_reco",&vtxx_reco,"vtxx_reco/D");
  fTree->Branch("vtxy_reco",&vtxy_reco,"vtxy_reco/D");
  fTree->Branch("vtxz_reco",&vtxz_reco,"vtxz_reco/D");
  fTree->Branch("trackstart_dcosx_reco",&trackstart_dcosx_reco, "trackstart_dcosx_reco/D");//
  fTree->Branch("trackstart_dcosy_reco",&trackstart_dcosy_reco, "trackstart_dcosy_reco/D");//
  fTree->Branch("trackstart_dcosz_reco",&trackstart_dcosz_reco, "trackstart_dcosz_reco/D");//
  fTree->Branch("trackexit_dcosx_reco",&trackexit_dcosx_reco, "trackexit_dcosx_reco/D");//
  fTree->Branch("trackexit_dcosy_reco",&trackexit_dcosy_reco, "trackexit_dcosy_reco/D");//
  fTree->Branch("trackexit_dcosz_reco",&trackexit_dcosz_reco, "trackexit_dcosz_reco/D");//
  fTree->Branch("trackstart_x_reco",&trackstart_x_reco, "trackstart_x_reco/D");
  fTree->Branch("trackstart_y_reco",&trackstart_y_reco, "trackstart_y_reco/D");
  fTree->Branch("trackstart_z_reco",&trackstart_z_reco, "trackstart_z_reco/D");
  fTree->Branch("trackexit_x_reco",&trackexit_x_reco, "trackexit_x_reco/D");
  fTree->Branch("trackexit_y_reco",&trackexit_y_reco, "trackexit_y_reco/D");
  fTree->Branch("trackexit_z_reco",&trackexit_z_reco, "trackexit_z_reco/D");    
  fTree->Branch("nmatched_reco",&nmatched_reco,"nmatched_reco/I");  
  fTree->Branch("trk_mom_minos",&trk_mom_minos,"trk_mom_minos/D");
  fTree->Branch("trk_charge_minos",&trk_charge_minos,"trk_charge_minos/D");
  fTree->Branch("trk_dcosx_minos",&trk_dcosx_minos,"trk_dcosx_minos/D");
  fTree->Branch("trk_dcosy_minos",&trk_dcosy_minos,"trk_dcosy_minos/D");
  fTree->Branch("trk_dcosz_minos",&trk_dcosz_minos,"trk_dcosz_minos/D");
  fTree->Branch("trk_vtxx_minos",&trk_vtxx_minos,"trk_vtxx_minos/D");
  fTree->Branch("trk_vtxy_minos",&trk_vtxy_minos,"trk_vtxy_minos/D");
  fTree->Branch("trk_vtxz_minos",&trk_vtxz_minos,"trk_vtxz_minos/D"); 
  fTree->Branch("test_charge_minos",&test_charge_minos,"test_charge_minos/I"); 
  fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  fTree->Branch("enu_truth",&enu_truth,"enu_truth/D");
  fTree->Branch("Q2_truth",&Q2_truth,"Q2_truth/D");
  fTree->Branch("hitnuc_truth",&hitnuc_truth,"hitnuc_truth/I");
  fTree->Branch("W_truth",&W_truth,"W_truth/D");
  fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/D");
  fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/D");
  fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/D");
  fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
  fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
  fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
  fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");

  fTree->Branch("nminos_tracks",&nminos_tracks,"nminos_tracks/I");  
  fTree->Branch("trk_charge_minos_all",trk_charge_minos_all,"trk_charge_minos_all[nminos_tracks]/D");
  fTree->Branch("muon_reco",&muon_reco,"muon_reco/I"); 
  fTree->Branch("minos_enter_true",&minos_enter_true,"minos_enter_true/I");

  fTree->Branch("trackstart_x_reco_muon",&trackstart_x_reco_muon, "trackstart_x_reco_muon/D");
  fTree->Branch("trackstart_y_reco_muon",&trackstart_y_reco_muon, "trackstart_y_reco_muon/D");
  fTree->Branch("trackstart_z_reco_muon",&trackstart_z_reco_muon, "trackstart_z_reco_muon/D");
  fTree->Branch("trackexit_x_reco_muon",&trackexit_x_reco_muon, "trackexit_x_reco_muon/D");
  fTree->Branch("trackexit_y_reco_muon",&trackexit_y_reco_muon, "trackexit_y_reco_muon/D");
  fTree->Branch("trackexit_z_reco_muon",&trackexit_z_reco_muon, "trackexit_z_reco_muon/D");    
  
  fTree->Branch("muon_exits",&muon_exits,"muon_exits/I"); 
}


void t962::CCRecoEffTree::beginSubRun(const art::SubRun& sr)
{

art::Handle< sumdata::POTSummary > potListHandle;
sr.getByLabel(fPOTModuleLabel,potListHandle);

if(sr.getByLabel(fPOTModuleLabel,potListHandle))
pot=potListHandle->totpot;
else
pot=0.;

}

void t962::CCRecoEffTree::analyze(const art::Event& evt)
{

  art::ServiceHandle<cheat::BackTracker> bt;
  bt->SetEveIdCalculator(new sim::EmEveIdCalculator);

  ResetVars();

  run = evt.run();
  event = evt.id().event();

  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;

  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  art::Handle< std::vector<t962::MINOS> > minosListHandle;
  evt.getByLabel(fMINOSModuleLabel,minosListHandle);
  art::Handle< std::vector<sim::Particle> > larListHandle;
  evt.getByLabel (fLArG4ModuleLabel,larListHandle);  

  art::PtrVector<simb::MCTruth> mclist;
  if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    }

  art::PtrVector<recob::Track> tracklist;
  if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
  for (unsigned int i = 0; i < trackListHandle->size(); ++i){
    art::Ptr<recob::Track> trackHolder(trackListHandle,i);
    tracklist.push_back(trackHolder);
  }

  art::PtrVector<recob::Vertex> vertexlist;
  if(evt.getByLabel(fVertexModuleLabel,vertexListHandle))
  for (unsigned int i = 0; i < vertexListHandle->size(); ++i){
    art::Ptr<recob::Vertex> vertexHolder(vertexListHandle,i);
    vertexlist.push_back(vertexHolder);
  }

  art::PtrVector<t962::MINOS> minoslist;
  if(evt.getByLabel(fMINOSModuleLabel,minosListHandle))
  for (unsigned int i = 0; i < minosListHandle->size(); i++){
    art::Ptr<t962::MINOS> minosHolder(minosListHandle,i);
    minoslist.push_back(minosHolder);
  }
  
  art::PtrVector<sim::Particle> parIn;
  if(evt.getByLabel(fLArG4ModuleLabel,larListHandle))
  for(unsigned int i = 0; i < larListHandle->size(); ++i){
      art::Ptr<sim::Particle> particle(larListHandle, i);
      parIn.push_back(particle);
  }

  art::ServiceHandle<geo::Geometry> geom;

  int MuonID = 0;
  unsigned int traj_points_muon = 0;;
  simb::MCTrajectory traj_muon;
  
  for(unsigned int j = 0; j < parIn.size(); ++j){ 
    if(abs(parIn[j]->PdgCode())== 13 && parIn[j]->Process()=="primary"){
      MuonID = parIn[j]->TrackId();
      traj_muon = parIn[j]->Trajectory();
      traj_points_muon = parIn[j]->NumberTrajectoryPoints();
    }
  }
  
  if(traj_points_muon > 1){
    
    for(unsigned int i = 1; i < traj_points_muon; i++){
      double Z1 = traj_muon.Position(i-1).Z();
      double Z2 = traj_muon.Position(i).Z();
      
      if(Z1 < 152.95 && Z2 >= 152.95){
	minos_enter_true = 1;
      }
    }
  } // if traj_points_muon > 1
  
  
    // Now see the track_reco  hits
  
  std::vector<int> skip_next_mu;
  std::vector<double> trackStart_muon;
  std::vector<double> trackEnd_muon;
  trackStart_muon.clear();
  trackEnd_muon.clear();
  
  art::FindManyP<recob::Hit> fmh(trackListHandle, evt, fTrackModuleLabel);

  for(unsigned int i = 0; i<tracklist.size(); ++i){
    
    std::vector< art::Ptr<recob::Hit> > track_reco_hits;
    double muon_hits_reco_track = 0;

    std::vector< art::Ptr<recob::Hit> > track_hitlist = fmh.at(i);
        
    track_reco_hits = track_hitlist;
    
    // loop over the hits and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Hit> >::iterator track_reco_itr = track_reco_hits.begin();
    
    
    while( track_reco_itr != track_reco_hits.end() ){

      std::vector<cheat::TrackIDE> track_reco_eveides   = bt->HitToEveID(*track_reco_itr);
      std::vector<cheat::TrackIDE>::iterator track_reco_idesitr = track_reco_eveides.begin();
      while( track_reco_idesitr != track_reco_eveides.end() ){
	if((*track_reco_idesitr).trackID == MuonID)
	  muon_hits_reco_track += (*track_reco_idesitr).energyFrac ;
	
	track_reco_idesitr++;	
      }      
      track_reco_itr++;     
    }
    
    if((muon_hits_reco_track/track_reco_hits.size())>=0.7){
      if(skip_next_mu.size()==0){
	skip_next_mu.push_back(1);
	muon_reco = 1;
	tracklist[i]->Extent(trackStart_muon,trackEnd_muon);
	trackstart_x_reco_muon=trackStart_muon[0];
	trackstart_y_reco_muon=trackStart_muon[1];
	trackstart_z_reco_muon=trackStart_muon[2];
	trackexit_x_reco_muon=trackEnd_muon[0];
	trackexit_y_reco_muon=trackEnd_muon[1];
	trackexit_z_reco_muon=trackEnd_muon[2];
	
	if(EndsOnBoundary(tracklist[i]))
	  muon_exits = 1; 
      }
    }
  }// loop over tracks
  
  //vertex information
  if(vertexlist.size())
    {
      double vtxxyz[3];
      vertexlist[0]->XYZ(vtxxyz);
      vtxx_reco = vtxxyz[0];
      vtxy_reco = vtxxyz[1];
      vtxz_reco = vtxxyz[2];
    }
  
  //matching information  

  //find matched MINOS information for each track
  art::FindOne<t962::MINOS> fomatch(trackListHandle, evt, fTrackMatchModuleLabel);

 
  test_charge_minos=0.;

  nminos_tracks = minoslist.size();
  
  for(unsigned int j = 0; j < minoslist.size(); j++){ 
    if (!isdata){
      trk_charge_minos_all[j] = minoslist[j]->fcharge;
      if(minoslist[j]->fcharge<0.)
	test_charge_minos=-1;
    }
  }
  
 
  for(unsigned int i = 0; i < tracklist.size(); i++)
  {
     if(!fomatch.at(i).isValid()) continue;//No matching MINOS track
     
     if(fomatch.at(i).ref().ftrkcontained)
        trk_mom_minos = fomatch.at(i).ref().ftrkErange;
     else
        trk_mom_minos = fomatch.at(i).ref().ftrkmom;     
     
     trk_charge_minos = fomatch.at(i).ref().fcharge;
     trk_dcosx_minos = fomatch.at(i).ref().ftrkdcosx;
     trk_dcosy_minos = fomatch.at(i).ref().ftrkdcosy;
     trk_dcosz_minos = fomatch.at(i).ref().ftrkdcosz;
     trk_vtxx_minos = fomatch.at(i).ref().ftrkVtxX;
     trk_vtxy_minos = fomatch.at(i).ref().ftrkVtxY;
     trk_vtxz_minos = fomatch.at(i).ref().ftrkVtxZ;     
  }    
      
  
  
  //track information
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  trackStart.clear();
  trackEnd.clear();
  
  for(unsigned int i=0; i<tracklist.size();++i){
    tracklist[i]->Extent(trackStart,trackEnd); 
    trackstart_x_reco=trackStart[0];
    trackstart_y_reco=trackStart[1];
    trackstart_z_reco=trackStart[2];
    trackexit_x_reco=trackEnd[0];
    trackexit_y_reco=trackEnd[1];
    trackexit_z_reco=trackEnd[2];   
  }
  
  trackStart.clear();
  trackEnd.clear();
  memset(larStart, 0, 3);
  memset(larEnd, 0, 3);

  for(unsigned int i=0; i<tracklist.size();++i){

    if(!fomatch.at(i).isValid()) continue;//No matching MINOS track
    ++nmatched_reco;
    tracklist[i]->Direction(larStart,larEnd);
    tracklist[i]->Extent(trackStart,trackEnd); 
     
    trackstart_dcosx_reco = larStart[0];
    trackstart_dcosy_reco = larStart[1];
    trackstart_dcosz_reco = larStart[2];       
    trackexit_dcosx_reco = larEnd[0];
    trackexit_dcosy_reco = larEnd[1];
    trackexit_dcosz_reco = larEnd[2];
    
    trackstart_x_reco=trackStart[0];
    trackstart_y_reco=trackStart[1];
    trackstart_z_reco=trackStart[2];
    trackexit_x_reco=trackEnd[0];
    trackexit_y_reco=trackEnd[1];
    trackexit_z_reco=trackEnd[2];  
  }
  

    //mc truth information
  if (!isdata){
    nuPDG_truth = mclist[0]->GetNeutrino().Nu().PdgCode();
    ccnc_truth = mclist[0]->GetNeutrino().CCNC();
    mode_truth = mclist[0]->GetNeutrino().Mode();
    Q2_truth = mclist[0]->GetNeutrino().QSqr();
    W_truth = mclist[0]->GetNeutrino().W();
    hitnuc_truth = mclist[0]->GetNeutrino().HitNuc();
    enu_truth = mclist[0]->GetNeutrino().Nu().E();
    nuvtxx_truth = mclist[0]->GetNeutrino().Nu().Vx();
    nuvtxy_truth = mclist[0]->GetNeutrino().Nu().Vy();
    nuvtxz_truth = mclist[0]->GetNeutrino().Nu().Vz();
    lep_mom_truth = mclist[0]->GetNeutrino().Lepton().P();
    if (mclist[0]->GetNeutrino().Lepton().P()){
      lep_dcosx_truth = mclist[0]->GetNeutrino().Lepton().Px()/mclist[0]->GetNeutrino().Lepton().P();
      lep_dcosy_truth = mclist[0]->GetNeutrino().Lepton().Py()/mclist[0]->GetNeutrino().Lepton().P();
      lep_dcosz_truth = mclist[0]->GetNeutrino().Lepton().Pz()/mclist[0]->GetNeutrino().Lepton().P();
    }
  }
  
  fTree->Fill();
}

  //---------------------------------------------------------------- 
void t962::CCRecoEffTree::ResetVars(){

  run = -99999;
  event = -99999;
  isdata = -99999;
  vtxx_reco = -99999;
  vtxy_reco  = -99999;
  vtxz_reco  = -99999;
  trackstart_x_reco = -99999;
  trackstart_y_reco = -99999;
  trackstart_z_reco = -99999;
  trackexit_x_reco = -99999;
  trackexit_y_reco = -99999;
  trackexit_z_reco = -99999;  
  trackstart_dcosx_reco = -99999;
  trackstart_dcosy_reco = -99999;
  trackstart_dcosz_reco = -99999;       
  trackexit_dcosx_reco = -99999;
  trackexit_dcosy_reco = -99999;
  trackexit_dcosz_reco = -99999;
  nmatched_reco = -99999;
  trk_mom_minos = -99999;
  trk_charge_minos = -99999;
  trk_dcosx_minos = -99999;
  trk_dcosy_minos = -99999;
  trk_dcosz_minos = -99999;  
  trk_vtxx_minos = -99999;
  trk_vtxy_minos = -99999;
  trk_vtxz_minos = -99999;
  
  test_charge_minos=-99999;

  nuPDG_truth = -99999;
  ccnc_truth = -99999;
  mode_truth = -99999;
  enu_truth = -99999;
  Q2_truth = -99999;
  W_truth = -99999;
  hitnuc_truth = -99999;
  nuvtxx_truth = -99999;
  nuvtxy_truth = -99999;
  nuvtxz_truth = -99999;
  lep_mom_truth = -99999;
  lep_dcosx_truth = -99999;
  lep_dcosy_truth = -99999;
  lep_dcosz_truth = -99999;

  nminos_tracks = -99999;
  muon_reco = -99999;
  minos_enter_true = -99999;

  trackstart_x_reco_muon = -99999;
  trackstart_y_reco_muon = -99999;
  trackstart_z_reco_muon = -99999;
  trackexit_x_reco_muon = -99999;
  trackexit_y_reco_muon = -99999;
  trackexit_z_reco_muon = -99999;

  muon_exits = -99999;  
}

bool t962::CCRecoEffTree::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
{
      std::vector<double> larStart, larEnd;
      lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)

 	if(fabs(larEnd[0])<fboundaryWindow 
	|| fabs(47.-larEnd[0])<fboundaryWindow 
	|| fabs(larEnd[1]+20.)<fboundaryWindow
	|| fabs(20.-larEnd[1])<fboundaryWindow 
	|| fabs(larEnd[2])<fboundaryWindow 
	|| fabs(90.-larEnd[2])<fboundaryWindow  )   
	return true;  
      else return false;
}








