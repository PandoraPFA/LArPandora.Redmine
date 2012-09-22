////////////////////////////////////////////////////////////////////////
/// file LambdaAna.cxx
//
/// author saima@ksu.edu
// 
/// 
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"


#include "T962/LambdaAna/LambdaAna.h"
#include "MCCheater/BackTracker.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "TDatabasePDG.h"
#include "RecoBase/recobase.h"
#include "Geometry/Geometry.h" 
#include "Utilities/AssociationUtil.h"


// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


//__________________________________________________________________
bool sort_pred2(const std::pair<art::Ptr<recob::Track>,double>& left, const std::pair<art::Ptr<recob::Track>,double>& right)
{
  return left.second < right.second;
}

//___________________________________________________________________
namespace cchyp {
  

  LambdaAna::LambdaAna(fhicl::ParameterSet const& pset): 
    fGenieModuleLabel(pset.get< std::string >("GenieModuleLabel")),
    fLArG4ModuleLabel(pset.get< std::string >("LArG4ModuleLabel")),
    fHitModuleLabel(pset.get< std::string >("HitModuleLabel")),
    fDetSimModuleLabel(pset.get< std::string >("DetSimModuleLabel")),
    fLineMergerModuleLabel(pset.get< std::string >("LineMergerModuleLabel")),
    fTrackModuleLabel(pset.get< std::string >("TrackModuleLabel")),
    fVertexModuleLabel(pset.get< std::string >("VertexModuleLabel"))

  {
    
  }
  //____________________________________________________________________
    void LambdaAna::beginJob(){

     // get access to the TFile service
     art::ServiceHandle<art::TFileService> tfs;

     fXreco_Xtrue = tfs->make<TH1F>("fXreco_Xtrue", "Xreco - Xtrue", 100,-10,10);
     fYreco_Ytrue = tfs->make<TH1F>("fYreco_Ytrue", "Yreco - Ytrue", 100,-10,10);
     fZreco_Ztrue = tfs->make<TH1F>("fZreco_Ztrue", "Zreco - Ztrue", 100,-10,10);

     fXreco_Xmuon = tfs->make<TH1F>("fXreco_Xmuon", "Xreco - Xmuon", 100,-10,10);
     fYreco_Ymuon = tfs->make<TH1F>("fYreco_Ymuon", "Yreco - Ymuon", 100,-10,10);
     fZreco_Zmuon = tfs->make<TH1F>("fZreco_Zmuon", "Zreco - Zmuon", 100,-10,10);
     fMINDIST = tfs->make<TH1F>("fMINDIST", "Min Dist B/W 2 tracks", 100,0,10);

     fRtracksEnu = tfs->make<TH2F>("fRtracksEnu", ";Neutrino Energy; No of Reco Tracks", 500,0, 10, 500, 0, 10);
     fTtracksEnu = tfs->make<TH2F>("fTtracksEnu", ";Neutrino Energy; No of True Tracks", 500,0, 10, 500, 0, 10);
     fRtracksEp = tfs->make<TH2F>("fRtracksEp", ";Proton Energy; No of Reco Tracks", 500,0, 10, 500, 0, 10);
     fTtracksEp = tfs->make<TH2F>("fTtracksEp", ";Proton Energy; No of True Tracks", 500,0, 10, 500, 0, 10);
     fTRtracks =  tfs->make<TH2F>("fTRtracks", ";# of True Tracks; # of Reco Tracks", 500,0, 10, 500, 0, 10);


     art::ServiceHandle<geo::Geometry> geo;
     
     double x = 2.1*geo->DetHalfWidth();
     double y = 2.1*geo->DetHalfHeight();
     double z = 2.*geo->DetLength();
     int xdiv = TMath::Nint(2*x/5.);  
     int ydiv = TMath::Nint(2*y/5.);  
     int zdiv = TMath::Nint(2*z/5.);  

     std::cout << "******* x y z div = ***********" << xdiv << " " << ydiv << " " << zdiv << std::endl; 
     
     fRVertexX = tfs->make<TH1F>("fRVertexX", ";x (cm)", xdiv,  -x, x);
     fRVertexY = tfs->make<TH1F>("fRVertexY", ";y (cm)", ydiv,  -y, y);
     fRVertexZ = tfs->make<TH1F>("fRVertexZ", ";z (cm)", zdiv, -0.2*z, z);
     
     fRVertexXY = tfs->make<TH2F>("fRVertexXY", ";x (cm);y (cm)", xdiv,     -x, x, ydiv, -y, y);
     fRVertexXZ = tfs->make<TH2F>("fRVertexXZ", ";z (cm);x (cm)", zdiv, -0.2*z, z, xdiv, -x, x);
     fRVertexYZ = tfs->make<TH2F>("fRVertexYZ", ";z (cm);y (cm)", zdiv, -0.2*z, z, ydiv, -y, y);

     fRmuonX_TmuonX = tfs->make<TH1F>("fRmuonX_TmuonX", "Xreco_muon - Xtrue_muon", 100,-25,25);
     fRmuonY_TmuonY = tfs->make<TH1F>("fRmuonY_TmuonY", "Yreco_muon - Ytrue_muon", 100,-25,25);
     fRmuonZ_TmuonZ = tfs->make<TH1F>("fRmuonZ_TmuonZ", "Zreco_muon - Ztrue_muon", 100,-50,50);
     
     fTVertexX = tfs->make<TH1F>("fTVertexX", ";x (cm)", xdiv,  -x, x);
     fTVertexY = tfs->make<TH1F>("fTVertexY", ";y (cm)", ydiv,  -y, y);
     fTVertexZ = tfs->make<TH1F>("fTVertexZ", ";z (cm)", zdiv, -0.2*z, z);

     fXtrue_vs_Xreco =  tfs->make<TH2F>("fXtrue_vs_Xreco", ";X_true (cm);X_reco (cm)", 500,     0, x, 500, 0, x);
     fYtrue_vs_Yreco =  tfs->make<TH2F>("fYtrue_vs_Yreco", ";Y_true (cm);Y_reco (cm)", 500,     -y, y, 500, -y, y);
     fZtrue_vs_Zreco =  tfs->make<TH2F>("fZtrue_vs_Zreco", ";Z_true (cm);Z_reco (cm)", 500, -0.2*z, z, 500, -0.2*z, z);

     fEventsWith_greaterTHAN_2cmX = tfs->make<TH1F>("fEventsWith_greaterTHAN_2cmX", ";Event Number", 10000,  0, 1000);
     fEventsWith_greaterTHAN_2cmY = tfs->make<TH1F>("fEventsWith_greaterTHAN_2cmY", ";Event Number", 10000,  0, 1000);
     fEventsWith_greaterTHAN_2cmZ = tfs->make<TH1F>("fEventsWith_greaterTHAN_2cmZ", ";Event Number", 10000,  0, 1000);

     fXtrue_vs_Xmuon = tfs->make<TH2F>("fXtrue_vs_Xmuon", ";X_true (cm);X_muon (cm)", 500,     0, x, 500, 0, x);
     fYtrue_vs_Ymuon = tfs->make<TH2F>("fYtrue_vs_Ymuon", ";Y_true (cm);Y_muon (cm)", 500,     -y, y, 500, -y, y);
     fZtrue_vs_Zmuon = tfs->make<TH2F>("fZtrue_vs_Zmuon", ";Z_true (cm);Z_muon (cm)", 500, -0.2*z, z, 500, -0.2*z, z);
     
     fXmuon_vs_Xreco = tfs->make<TH2F>("fXmuon_vs_Xreco", ";X_muon (cm);X_reco (cm)", 500,     0, x, 500, 0, x);
     fYmuon_vs_Yreco = tfs->make<TH2F>("fYmuon_vs_Yreco", ";Y_muon (cm);Y_reco (cm)", 500,     -y, y, 500, -y, y);
     fZmuon_vs_Zreco = tfs->make<TH2F>("fZmuon_vs_Zreco", ";Z_muon (cm);Z_reco (cm)", 500, -0.2*z, z, 500, -0.2*z, z);

     fnumber_true =  tfs->make<TH1F>("fnumbertrue", ";x (cm)", xdiv,  -x, x);

     fTotal_CCQE =  tfs->make<TH1F>("fTotal_CCQE", ";CCQE Hyperon Total", 500,  0, 1000);
     fCCQE_FV =  tfs->make<TH1F>("fCCQE_FV", ";CCQE Hyperon after FV cut on Primary", 500,  0, 1000);
     fCCQE_decay_FV =  tfs->make<TH1F>("fCCQE_decay_FV", ";CCQE Hyperon after FV cut on Primary AND secondary", 500,  0, 1000);
     fMuon_Escapes =  tfs->make<TH1F>("fMuon_Escapes", ";Muon_Escapes", 0.5,  0, 3);
     fMuon_Enters_MINOS =  tfs->make<TH1F>("fMuon_Enters_MINOS", ";Muon_Enters_MINOS", 0.5,  0, 3);
     fProton_Escapes =  tfs->make<TH1F>("fProton_Escapes", ";Proton_Escapes", 0.5,  0, 3);
     fPion_Escapes =  tfs->make<TH1F>("fPion_Escapes", ";Pion_Escapes", 0.5,  0, 3);

     fMuon_Length =  tfs->make<TH1F>("fMuon_Length", ";Muon True Length in TPC (cm)", 500,  0, 115); 
     fProton_Length = tfs->make<TH1F>("fProtonn_Length", ";Proton True Length in TPC (cm)", 500,  0, 115); 
     fPion_Length = tfs->make<TH1F>("fPion_Length", ";Pion True Length in TPC (cm)", 500,  0, 115);

     fMuon_KE_Length = tfs->make<TH2F>("fMuon_KE_Length", ";Muon KE (GeV);Muon True Length in TPC (cm)", 500,     0, 15, 500, 0, 115);
     fProton_KE_Length = tfs->make<TH2F>("fProton_KE_Length", ";Proton KE (GeV);Proton True Length in TPC (cm)", 500,     0, 3, 500, 0, 115);
     fPion_KE_Length = tfs->make<TH2F>("fPion_KE_Length", ";Pion KE (GeV);Pion True Length in TPC (cm)", 500,     0, 2, 500, 0, 115);
    
     fU_no_clus =  tfs->make<TH1F>("fU_no_clus", ";Number of Clusters in Induction", 100,  0, 10);
     fV_no_clus =  tfs->make<TH1F>("fV_no_clus", ";Number of Clusters in Collection", 100,  0, 10);

     fU_clu_length = tfs->make<TH1F>("fU_clu_length", ";Cluster Length in Induction", 100,  0, 1000);
     fV_clu_length = tfs->make<TH1F>("fV_clu_length", ";Cluster Length in Collection", 100,  0, 1000);
     fUwire_span_length = tfs->make<TH1F>("fUwire_span_length", ";Projection on Induction Plane (cm)", 100,  0, 100);
     fVwire_span_length = tfs->make<TH1F>("fVwire_span_length", ";Projection on Collection Plane (cm)", 100,  0, 100);

     fMuon_Reco_Eff = tfs->make<TH2F>("fMuon_Reco_Eff", ";Muon KE (GeV);Efficiency", 500,     0, 10, 500, 0, 1.5);
     fProtonn_Reco_Eff = tfs->make<TH2F>("fProton_Reco_Eff", ";Proton KE (GeV); Efficiency", 500,     0, 3, 500, 0, 1.5);
     fPion_Reco_Eff = tfs->make<TH2F>("fPion_Reco_Eff", ";Pion KE (GeV);Efficiency", 500,     0, 2, 500, 0, 1.5);

     fMuon_Reco_Eff_vs_hits = tfs->make<TH2F>("fMuon_Reco_Eff_vs_hits", ";Muon Cheated Hits;Efficiency", 500,     0, 500, 500, 0, 1.5);
     fProtonn_Reco_Eff_vs_hits = tfs->make<TH2F>("fProton_Reco_Eff_vs_hits", ";Proton Cheated Hits; Efficiency", 500,     0, 500, 500, 0, 1.5);
     fPion_Reco_Eff_vs_hits = tfs->make<TH2F>("fPion_Reco_Eff_vs_hits", ";Pion Cheated Hits;Efficiency", 500,     0, 500, 500, 0, 1.5);

     flambda_decay_length =  tfs->make<TH1F>("flambda_decay_length", ";Lambda Decay Length (cm)", 10,  0, 10);
     fmuon_lambda_angle = tfs->make<TH1F>("fmuon_lambda_angle", ";Muon_Lambda cos(angle)", 10,  0, 1);
     fproton_pion_angle =  tfs->make<TH1F>("fproton_pion_angle", ";Proton_Pion cos(angle)", 10,  0, 1);
     flambda_KE = tfs->make<TH1F>("flambda_KE", ";Lambda KE (GeV)", 10,  0, 3);
     
     fMuon_Reco_Eff_3d = tfs->make<TH2F>("fMuon_Reco_Eff_3d", ";Muon KE (GeV);Efficiency", 500,     0, 10, 500, 0, 1.5);
     fProtonn_Reco_Eff_3d = tfs->make<TH2F>("fProton_Reco_Eff_3d", ";Proton KE (GeV);Efficiency", 500,     0, 3, 500, 0, 1.5);
     fPion_Reco_Eff_3d = tfs->make<TH2F>("fPion_Reco_Eff_3d", ";Pion KE (GeV);Efficiency", 500,     0, 2, 500, 0, 1.5);
     fNumber_tracks = tfs->make<TH1F>("fNumber_tracks", ";Number of 3D tracks", 10,  0, 5);

     fpurity_muon = tfs->make<TH2F>("fMuon_Reco_Purity", ";Muon KE (GeV);Purity", 500,     0, 10, 500, 0, 1.5);
     fpurity_proton = tfs->make<TH2F>("fProton_Reco_Purity", ";Proton KE (GeV);Purity", 500,     0, 3, 500, 0, 1.5);
     fpurity_pion = tfs->make<TH2F>("fPion_Reco_Purity", ";Pion KE (GeV);Purity", 500,     0, 3, 500, 0, 1.5);

     fcompleteness_muon =  tfs->make<TH2F>("fMuon_Reco_Completeness", ";Muon KE (GeV);Completeness", 500,     0, 10, 500, 0, 1.5);
     fcompleteness_proton = tfs->make<TH2F>("fProton_Reco_Completeness", ";Proton KE (GeV);Completeness", 500,     0, 3, 500, 0, 1.5);
     fcompleteness_pion = tfs->make<TH2F>("fPion_Reco_Completeness", ";Pion KE (GeV);Completeness", 500,     0, 3, 500, 0, 1.5);

    fmuon_cluster_hits_UV =  tfs->make<TH2F>("fmuon_cluster_hits_UV", ";Muon Cluster Hits (Induction);Muon Cluster Hits (Collection)", 500,0,250,500,0,250);
    fproton_cluster_hits_UV =  tfs->make<TH2F>("fproton_cluster_hits_UV", ";Proton Cluster Hits (Induction);Proton Cluster Hits (Collection)", 500,0,150,500,0,150);
    fpion_cluster_hits_UV =  tfs->make<TH2F>("fpion_cluster_hits_UV", ";Pion Cluster Hits (Induction);Pion Cluster Hits (Collection)", 500,0,100,500,0,100);

    fmuon_reco_percentage =  tfs->make<TH1F>("fmuon_reco_percentage", ";muon reco: 0 = Not Reconstructed, 1 = Reconstructed", 7,  -1, 2.5);
    fproton_reco_percentage =  tfs->make<TH1F>("fproton_reco_percentage", ";proton reco: 0 = Not Reconstructed, 1 = Reconstructed", 7,  -1, 2.5);
    fpion_reco_percentage =  tfs->make<TH1F>("fpion_reco_percentage", ";pion reco: 0 = Not Reconstructed, 1 = Reconstructed", 7,  -1, 2.5);

    fmuon_cluster_charge_UV=tfs->make<TH2F>("fmuon_cluster_charge_UV", ";Muon Cluster Charge x 10^3 (Induction);Muon Cluster Charge 10^3 (Collection)", 500,0,200,500,0,200);
    fproton_cluster_charge_UV=tfs->make<TH2F>("fproton_cluster_charge_UV", ";Proton Cluster Charge x 10^3 (Induction);Proton Cluster Charge x 10^3(Collection)", 500,0,200,500,0,200);
    fpion_cluster_charge_UV=tfs->make<TH2F>("fpion_cluster_charge_UV", ";Pion Cluster Charge x 10^3 (Induction);Pion Cluster Charge x 10^3 (Collection)", 500,0,200,500,0,200);

  }

  //___________________________________________________________________
  LambdaAna::~LambdaAna()
  {

  }

  //__________________________________________________________________
  void LambdaAna::analyze(const art::Event& evt)
  {
    art::ServiceHandle<cheat::BackTracker> bt;
    
    mf::LogVerbatim("LamdaAna")<<"---------------------------------------------------"
			       <<"---------------------------\n"
			       << "event  : " << evt.id().event();

    
    // get the truth info 

    TLorentzVector vertex;
    //double energynu = 0.;

    art::Handle< std::vector<simb::MCTruth> > mclist;
    evt.getByLabel(fGenieModuleLabel,mclist);
    art::Ptr<simb::MCTruth> mc(mclist, 0);
    const simb::MCParticle nu = mc->GetNeutrino().Nu();
    vertex  = nu.Position(); ///////////////////////////// REQ INFO
//     energynu  = nu.E(); ////////////////////////////////// REQ INFO

   
//     const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
 
//     int ntracks = 0;
//     double protonE = 0.;
//     int nprotons = 0;
//     double singleproton = 0;

//     ///Loop over the particle stack for this event 
//     for(int i = 0; i < mc->NParticles(); ++i){
//       simb::MCParticle part(mc->GetParticle(i));

//       if(part.StatusCode() == 1 && (part.PdgCode() == 2212 || part.PdgCode() == 13)){ //finalstableparticle proton or muon
// 	ntracks += 1; ////////////////////////////////// REQ INFO
// 	//std::cout << "final state particle energy = " << part.E() << std::endl;
//       }

//       if(part.StatusCode() == 14 && part.PdgCode() == 2212){ // proton/hadron in the nucleus
// 	protonE = part.E(); ////////////////////////////////// REQ INFO
// 	//std::cout << "hadron in the nucleaus has energy = " << part.E() << std::endl;
//       }

//       if(part.StatusCode() == 1 && (part.PdgCode() == 2212)){ //finalstableparticle proton
// 	nprotons += 1;
//       }

//       if(nprotons==1){
// 	for(int i = 0; i < mc->NParticles(); ++i){
// 	  if(part.StatusCode() == 1 && (part.PdgCode() == 2212)){ //finalstableparticle proton
// 	    singleproton = part.E();
// 	    //std::cout << "SINGLE PROTON EVENT, with proton energy = " << singleproton << std::endl;
// 	    break;
// 	  }
// 	}
//       }

//     }
     std::cout << "the true vertex position = " << vertex.X() << " " << vertex.Y() << " " << vertex.Z() << std::endl;
//     std::cout << "number of true tracks in this event = " << ntracks << std::endl;
//     std::cout << "neutrino energy = " << energynu << std::endl;
//     std::cout << "proton energy = " << protonE << std::endl;



//------------------------------------

    //get the reco track info
    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrackModuleLabel,trackListHandle);
    
    art::PtrVector<recob::Track> trkIn;
    for(size_t ii = 0; ii < trackListHandle->size(); ++ii){
      art::Ptr<recob::Track> track(trackListHandle, ii);
      trkIn.push_back(track);
    }

    art::FindManyP<recob::Hit> fmht(trackListHandle, evt, fTrackModuleLabel);

//     std::vector<double> start;
//     std::vector<double> end;

//     TVector3 startXYZ;
//     TVector3 endXYZ;

//     std::vector< std::pair<art::Ptr<recob::Track>, double> > trackpair;
    
//     for(unsigned int i = 0; i<trkIn.size(); ++i){
//       trkIn[i]->Extent(start, end);
//       startXYZ.SetXYZ(start[0],start[1],start[2]);
//       endXYZ.SetXYZ(end[0],end[1],end[2]);


//       double length = (endXYZ-startXYZ).Mag();
//       //std::cout << "Track length calculated = " << length << std::endl;
//       trackpair.push_back(std::pair<art::Ptr<recob::Track>,double>(trkIn[i],length));
//     }
    
//     for(unsigned int i = 0; i<trackpair.size(); ++i){
//       std::cout << "track id is  = " << (trackpair[i].first)->ID() << " track length = " << (trackpair[i].second) << std::endl;
//     }
    
//     std::sort(trackpair.rbegin(), trackpair.rend(),sort_pred2);
    
//     std::cout << "AFTER SORTING " << std::endl;
//     for(unsigned int i = 0; i<trackpair.size(); ++i){
//       std::cout << "track id is  = " << (trackpair[i].first)->ID() << " track length = " << (trackpair[i].second) << std::endl;
//     }

//     std::vector<double> startpoint, endpoint;
//     TVector3 muon_start;   ////////////////////////////////////////// REQ INFO
//     int recotracks = trackpair.size();  ///////////////////////////// REQ INFO
    
//     if(recotracks>0){
//       trackpair[0].first->Extent(startpoint, endpoint);
//       muon_start.SetXYZ(startpoint[0],startpoint[1],startpoint[2]);
//     }

//     std::cout << "muon reco start = " << muon_start[0] << " " << muon_start[1] << " " << muon_start[2] << std::endl;
//     std::cout << "# of rec tracks = " << trkIn.size() << std::endl;

//---------------------------------------------

    //get the reco vertex info
    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    
    art::PtrVector<recob::Vertex> vtxIn;
    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
      {
	art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
	vtxIn.push_back(vertex);
      }

    double reco_X = 0.;
    double reco_Y = 0.;
    double reco_Z = 0.;
    Double_t vtxcoord[3];
    if (vtxIn.size() > 0){
      vtxIn[0]->XYZ(vtxcoord);

      reco_X = vtxcoord[0];//////////////////////// REQ INFO
      reco_Y = vtxcoord[1];//////////////////////// REQ INFO
      reco_Z = vtxcoord[2];//////////////////////// REQ INFO

    }

    std::cout << "reco vertex position = " << reco_X << " " << reco_Y << " " << reco_Z << std::endl;


    //// make histos now

//     fRtracksEnu->Fill(energynu, trkIn.size());
//     fTtracksEnu->Fill(energynu, ntracks);
//     fRtracksEp->Fill(protonE, trkIn.size());
//     fTtracksEp->Fill(protonE, ntracks);
//     fTRtracks->Fill(ntracks, trkIn.size());
    
    
//     if(vtxIn.size()>0 &&  
//        reco_X > 3 && reco_X < 44 &&
//        reco_Y > -16 && reco_Y < 16 &&
//        reco_Z > 6 && reco_Z < 86){
      if(vtxIn.size()>0){
      
// //       fXreco_Xtrue->Fill(reco_X-vertex.X());
// //       fYreco_Ytrue->Fill(reco_Y-vertex.Y()); 
// //       fZreco_Ztrue->Fill(reco_Z-vertex.Z());
      
// //       fXreco_Xmuon->Fill(reco_X-muon_start[0]);
// //       fYreco_Ymuon->Fill(reco_Y-muon_start[1]); 
// //       fZreco_Zmuon->Fill(reco_Z-muon_start[2]);
      
// //       fRVertexX->Fill(reco_X);
// //       fRVertexY->Fill(reco_Y);
// //       fRVertexZ->Fill(reco_Z);
      
// //       fRVertexXY->Fill(reco_X, reco_Y); //looking at the beam
// //       fRVertexXZ->Fill(reco_Z, reco_X); //top view
// //       fRVertexYZ->Fill(reco_Z, reco_Y); //side view
      
// //       fRmuonX_TmuonX->Fill(muon_start[0]-vertex.X());
// //       fRmuonY_TmuonY->Fill(muon_start[1]-vertex.Y());	 
// //       fRmuonZ_TmuonZ->Fill(muon_start[2]-vertex.Z());
      
// //       fTVertexX->Fill(vertex.X());
// //       fTVertexY->Fill(vertex.Y());
// //       fTVertexZ->Fill(vertex.Z());

      fXtrue_vs_Xreco->Fill(vertex.X(), reco_X);
      fYtrue_vs_Yreco->Fill(vertex.Y(), reco_Y);
      fZtrue_vs_Zreco->Fill(vertex.Z(), reco_Z);
 
//       if(abs(reco_X-vertex.X()) > 2){
// 	fEventsWith_greaterTHAN_2cmX->Fill(evt.id().event());
// 	std::cout << "for reco_X-vertex.X() > 2 , events are " << evt.id().event() << std:: endl;
//       }

//       if(abs(reco_Y-vertex.Y()) > 2){
// 	fEventsWith_greaterTHAN_2cmY->Fill(evt.id().event()); 
//       	std::cout << "for reco_Y-vertex.Y() > 2 , events are " << evt.id().event() << std:: endl;
//       }

//       if(abs(reco_Z-vertex.Z()) > 2){
// 	fEventsWith_greaterTHAN_2cmZ->Fill(evt.id().event()); 
// 	std::cout << "for reco_Z-vertex.Z() > 2 , events are " << evt.id().event() << std:: endl;
//       }

//       fXtrue_vs_Xmuon->Fill(vertex.X(), muon_start[0]);
//       fYtrue_vs_Ymuon->Fill(vertex.Y(), muon_start[1]);
//       fZtrue_vs_Zmuon->Fill(vertex.Z(), muon_start[2]);

//       fXmuon_vs_Xreco->Fill(muon_start[0], reco_X);
//       fYmuon_vs_Yreco->Fill(muon_start[1], reco_Y);
//       fZmuon_vs_Zreco->Fill(muon_start[2], reco_Z);

//       if(vertex.X() > 3 && vertex.X() < 44 &&
// 	 vertex.Y() > -16 && vertex.Y() < 16 &&
// 	 vertex.Z() > 6 && vertex.Z() < 86 ) {

// 	fnumber_true->Fill(vertex.X());
//       }

      
      }
//______________________________________________
      TLorentzVector lMuonStart;
      TLorentzVector lProtonStart;
      TLorentzVector lPionStart;

      TLorentzVector lMuonMom;
      TLorentzVector lProtonMom;
      TLorentzVector lPionMom;

      TLorentzVector lLambdaMom;
      TVector3 LambdaMom;
      double LambdaMass = 0.;
      
      TVector3 MuonStart;
      TVector3 ProtonStart;
      TVector3 PionStart;

      TVector3 MuonEnd;
      TVector3 ProtonEnd;
      TVector3 PionEnd;

      TVector3 MuonMom;
      TVector3 ProtonMom;
      TVector3 PionMom;

      //double MuonE = 0.;
      //double ProtonE = 0.;
      //double PionE = 0.;

      int MuonID = 0;
      int ProtonID = 0;
      int PionID = 0;

      double MuonMass = 0.;
      double ProtonMass = 0.;
      double PionMass = 0.;

      double MuonLength = 0.;
      double ProtonLength = 0.;
      double PionLength = 0.;

      double Muon_theta = 0.;
      double Muon_phi = 0.;
      double Proton_theta = 0;
      double Proton_phi = 0;
      double Pion_theta = 0;
      double Pion_phi = 0;

      simb::MCTrajectory traj_muon;
      unsigned int traj_points_muon; 

      simb::MCTrajectory traj_proton;
      unsigned int traj_points_proton; 

      simb::MCTrajectory traj_pion;
      unsigned int traj_points_pion; 

      
    art::Handle< std::vector<sim::Particle> > larListHandle;
    evt.getByLabel (fLArG4ModuleLabel,larListHandle);  

    art::PtrVector<sim::Particle> parIn;
    for(size_t i = 0; i < larListHandle->size(); ++i){
      art::Ptr<sim::Particle> particle(larListHandle, i);
      parIn.push_back(particle);
    }


    for(size_t a = 0; a < parIn.size(); ++a){
      if(abs(parIn[a]->PdgCode())==13 && parIn[a]->Process()=="primary"){  // if CCQE numubar event
	
	for(size_t b = 0; b < parIn.size(); ++b){
	  if(abs(parIn[b]->PdgCode())==3122 && parIn[b]->Process()=="primary"){ 
	    // if Lambda is in the Primary Particles
	    int LID = parIn[b]->TrackId();

	    for(unsigned int c = 0; c < parIn.size(); ++c){ 
	      if(parIn[c]->Mother()== LID && parIn[c]->Process()=="Decay" && 
		 (parIn[c]->PdgCode()==2212)){
		// if lambda decays to proton
		//break;
		int EventNumber =  evt.id().event();
		mf::LogInfo("LambdaAna") << "event# = " << EventNumber;
		
		fTotal_CCQE->Fill(EventNumber);
		
		
		for(unsigned int d = 0; d < parIn.size(); ++d){ // if primary vertex in FV
		  if(abs(parIn[d]->PdgCode())==13 && parIn[d]->Process()=="primary" 
		     && parIn[d]->Vx()>0   && parIn[d]->Vx()<47 
		     && parIn[d]->Vy()>-20 && parIn[d]->Vy()<20 
		     && parIn[d]->Vz()>0   && parIn[d]->Vz()<90){ 
		    
		    fCCQE_FV->Fill(evt.id().event());
		    
		    for(unsigned int e = 0; e < parIn.size(); ++e){ // if secondary vertex in FV
		      if(parIn[e]->Mother()== LID && parIn[e]->Process()=="Decay" && (parIn[e]->PdgCode()==2212)
			 && parIn[e]->Vx()>0   && parIn[d]->Vx()<47 
			 && parIn[e]->Vy()>-20 && parIn[d]->Vy()<20 
			 && parIn[e]->Vz()>0   && parIn[d]->Vz()<90){ 
			
			fCCQE_decay_FV->Fill(evt.id().event());
			
			// get mu+ info
			for(unsigned int j = 0; j < parIn.size(); ++j){ 
			  if(parIn[j]->PdgCode()== -13 && parIn[j]->Process()=="primary"){
			    
			    lMuonStart = parIn[j]->Position();
			    MuonStart = lMuonStart.Vect();
			    MuonEnd = parIn[j]->EndPoint();
			    lMuonMom = parIn[j]->Momentum();
			    MuonMom = lMuonMom.Vect();
			    MuonID = parIn[j]->TrackId();
			    MuonMass = parIn[j]->Mass();
			    traj_muon =  parIn[j]->Trajectory();
			    traj_points_muon = parIn[j]->NumberTrajectoryPoints(); 
			  }
			}
			
			// get the proton info
			for(unsigned int j = 0; j < parIn.size(); ++j){
			  if(parIn[j]->Mother()== LID && parIn[j]->Process()=="Decay" && (parIn[j]->PdgCode()==2212)){
			    
			    lProtonStart = parIn[j]->Position();
			    ProtonStart = lProtonStart.Vect();
			    ProtonEnd = parIn[j]->EndPoint();
			    lProtonMom = parIn[j]->Momentum();
			    ProtonMom = lProtonMom.Vect();
			    ProtonID = parIn[j]->TrackId();
			    ProtonMass = parIn[j]->Mass();
			    traj_proton =  parIn[j]->Trajectory();
			    traj_points_proton = parIn[j]->NumberTrajectoryPoints(); 
			  }
			}
			
			// get the pion info
			for(unsigned int j = 0; j < parIn.size(); ++j){
			  if(parIn[j]->Mother()== LID && parIn[j]->Process()=="Decay" && (parIn[j]->PdgCode()==-211)){
			    
			    lPionStart = parIn[j]->Position();
			    PionStart = lPionStart.Vect();
			    PionEnd = parIn[j]->EndPoint();
			    lPionMom = parIn[j]->Momentum();
			    PionMom = lPionMom.Vect();
			    PionID = parIn[j]->TrackId();
			    PionMass = parIn[j]->Mass();
			    traj_pion =  parIn[j]->Trajectory();
			    traj_points_pion = parIn[j]->NumberTrajectoryPoints(); 
			  }
			}
			
			// get lambda info - a little
			for(unsigned int f = 0; f < parIn.size(); ++f){
			  if(abs(parIn[f]->PdgCode())==3122 && parIn[f]->Process()=="primary"){
			    lLambdaMom = parIn[f]->Momentum();
			    LambdaMom = lLambdaMom.Vect();
			    LambdaMass = parIn[f]->Mass();
			  }
			}
			
			
			//find if the muon is exiting
			mf::LogInfo("LamdaAna") << "number of trajectory points of muon = " 
						<< traj_points_muon;
			
			for(unsigned int i = 1; i < traj_points_muon; ++i){
			  double X1 = traj_muon.X(i-1);
			  double Y1 = traj_muon.Y(i-1);
			  double Z1 = traj_muon.Z(i-1);
			  double X2 = traj_muon.X(i);
			  double Y2 = traj_muon.Y(i);
			  double Z2 = traj_muon.Z(i);
			  //std::cout << "POSITION of MUON along its way = " << X1 << " " << Y1 << " " << Z1 << std::endl;
			  if(  (((X1<47&&X2>=47) || (X1>0&&X2<=0)) && (Y1>=-20&&Y1<20) && (Z1>0&&Z1<=90)) 
			       || (((Y1<20&&Y2>=20) || (Y1>-20&&Y2<=-20)) && (X1>0&&X1<=47) && (Z1>0&&Z1<=90))
			       || (((Z1<90&&Z2>=90) || (Z1>0&&Z2<=0)) && (X1>0&&X1<=47) && (Y1>=-20&&Y1<20))
			       ){
			    
			    
			    TVector3 muon_TPC_exit(X2, Y2, Z2);
			    MuonEnd = muon_TPC_exit;
			    //std::cout << "*****************" << std::endl;
			    
			    fMuon_Escapes->Fill(1);
			  }
			  
			  if(Z1 < 152.95 && Z2 >= 152.95){
			    fMuon_Enters_MINOS->Fill(1);
			  }
			}
			
			//find if the proton is exiting
			mf::LogInfo("LamdaAna") << "number of trajectory points of proton = " 
						<< traj_points_proton;
			
			for(unsigned int i = 1; i < traj_points_proton; ++i){
			  double X1 = traj_proton.X(i-1);
			  double Y1 = traj_proton.Y(i-1);
			  double Z1 = traj_proton.Z(i-1);
			  double X2 = traj_proton.X(i);
			  double Y2 = traj_proton.Y(i);
			  double Z2 = traj_proton.Z(i);
			  //std::cout << "POSITION of PROTON along its way = " << X1 << " " << Y1 << " " << Z1 << std::endl;
			  if(  (((X1<47&&X2>=47) || (X1>0&&X2<=0)) && (Y1>=-20&&Y1<20) && (Z1>0&&Z1<=90)) 
			       || (((Y1<20&&Y2>=20) || (Y1>-20&&Y2<=-20)) && (X1>0&&X1<=47) && (Z1>0&&Z1<=90))
			       || (((Z1<90&&Z2>=90) || (Z1>0&&Z2<=0)) && (X1>0&&X1<=47) && (Y1>=-20&&Y1<20))
			       ){
			    
			    TVector3 proton_TPC_exit(X2, Y2, Z2);
			    ProtonEnd = proton_TPC_exit;
			    //std::cout << "*****************" << std::endl;
			    
			    fProton_Escapes->Fill(1);
			  }
			}
			
			//find if the pion is exiting
			std::cout << "number of trajectory points of pion = " << traj_points_pion << std::endl;
			
			for(unsigned int i = 1; i < traj_points_pion; ++i){
			  double X1 = traj_pion.X(i-1);
			  double Y1 = traj_pion.Y(i-1);
			  double Z1 = traj_pion.Z(i-1);
			  double X2 = traj_pion.X(i);
			  double Y2 = traj_pion.Y(i);
			  double Z2 = traj_pion.Z(i);
			  //std::cout << "POSITION of PION along its way = " << X1 << " " << Y1 << " " << Z1 << std::endl;
			  if(  (((X1<47&&X2>=47) || (X1>0&&X2<=0)) && (Y1>=-20&&Y1<20) && (Z1>0&&Z1<=90)) 
			       || (((Y1<20&&Y2>=20) || (Y1>-20&&Y2<=-20)) && (X1>0&&X1<=47) && (Z1>0&&Z1<=90))
			       || (((Z1<90&&Z2>=90) || (Z1>0&&Z2<=0)) && (X1>0&&X1<=47) && (Y1>=-20&&Y1<20))
			       ){
			    
			    TVector3 pion_TPC_exit(X2, Y2, Z2);
			    PionEnd = pion_TPC_exit;
			    //std::cout << "*****************" << std::endl;
			    
			    fPion_Escapes->Fill(1);
			  }
			}
			
			
			std::cout << "muon start, end and momentum = " << std::endl;
			MuonStart.Print();
			MuonEnd.Print();
			MuonMom.Print();
			std::cout << "muon ID = " << MuonID << " muon mass = " << MuonMass <<  std::endl;
			
			std::cout << "proton start, end and momentum = " << std::endl;
			ProtonStart.Print();
			ProtonEnd.Print();
			ProtonMom.Print();
			std::cout << "proton ID = " << ProtonID << " proton mass = " << ProtonMass<< std::endl;
			
			std::cout << "pion start, end and momentum = " << std::endl;
			PionStart.Print();
			PionEnd.Print();
			PionMom.Print();
			std::cout << "pion ID = " << PionID << " pion mass = " << PionMass<< std::endl;
			
			
			// some calculations
			MuonLength = (MuonEnd-MuonStart).Mag();
			ProtonLength  = (ProtonEnd-ProtonStart).Mag();
			PionLength = (PionEnd-PionStart).Mag();
			
			std::cout << "muon, proton and pion length = " << MuonLength << " " << ProtonLength << " " << PionLength << std::endl;
			
			// theta, phi calculation might not be correct or might not match what we get from root, need to work on this !! 
			Muon_theta = (TMath::ACos(MuonMom.Z()/MuonMom.Mag()))*(180/TMath::Pi());
			Proton_theta = (TMath::ACos(ProtonMom.Z()/ProtonMom.Mag()))*(180/TMath::Pi());
			Pion_theta = (TMath::ACos(PionMom.Z()/PionMom.Mag()))*(180/TMath::Pi());
			
			std::cout << "theta of Muon, proton, pion = " << Muon_theta << " " 
				  << Proton_theta << " " << Pion_theta << std::endl;
			
			if(MuonMom.X()>0 && MuonMom.Z()>0)
			  Muon_phi = (TMath::ATan(MuonMom.X()/MuonMom.Z()))*(180/TMath::Pi());
			
			if((MuonMom.X()>0 && MuonMom.Z()<0) ||(MuonMom.X()<0 && MuonMom.Z()<0))
			  Muon_phi = (TMath::Pi()+(TMath::ATan(MuonMom.X()/MuonMom.Z())))*(180/TMath::Pi());
			
			if(MuonMom.X()<0 && MuonMom.Z()>0)
			  Muon_phi = ((2*(TMath::Pi()))+(TMath::ATan(MuonMom.X()/MuonMom.Z())))*(180/TMath::Pi());
			
			if(ProtonMom.X()>0 && ProtonMom.Z()>0)
			  Proton_phi = (TMath::ATan(ProtonMom.X()/ProtonMom.Z()))*(180/TMath::Pi());
			
			if((ProtonMom.X()>0 && ProtonMom.Z()<0) ||(ProtonMom.X()<0 && ProtonMom.Z()<0))
			  Proton_phi = (TMath::Pi()+(TMath::ATan(ProtonMom.X()/ProtonMom.Z())))*(180/TMath::Pi());
			
			if(ProtonMom.X()<0 && ProtonMom.Z()>0)
			  Proton_phi = ((2*(TMath::Pi()))+(TMath::ATan(ProtonMom.X()/ProtonMom.Z())))*(180/TMath::Pi());
			
			
			if(PionMom.X()>0 && PionMom.Z()>0)
			  Pion_phi = (TMath::ATan(PionMom.X()/PionMom.Z()))*(180/TMath::Pi());
			
			if((PionMom.X()>0 && PionMom.Z()<0) ||(PionMom.X()<0 && PionMom.Z()<0))
			  Pion_phi = (TMath::Pi()+(TMath::ATan(PionMom.X()/PionMom.Z())))*(180/TMath::Pi());
			
			if(PionMom.X()<0 && PionMom.Z()>0)
			  Pion_phi = ((2*(TMath::Pi()))+(TMath::ATan(PionMom.X()/PionMom.Z())))*(180/TMath::Pi());
			
			std::cout << "phi of Muon, proton, pion = " << Muon_phi << " " 
				  << Proton_phi << " " << Pion_phi << std::endl;
			
			//fill the histos
			fMuon_Length->Fill(MuonLength);
			fProton_Length->Fill(ProtonLength);
			fPion_Length->Fill(PionLength);
			
			fMuon_KE_Length->Fill((lMuonMom.E()-(MuonMass/1000)), MuonLength);
			fProton_KE_Length->Fill((lProtonMom.E()-(ProtonMass/1000)), ProtonLength);
			fPion_KE_Length->Fill((lPionMom.E()-(PionMass/1000)), PionLength);
			
			flambda_decay_length->Fill((ProtonStart-MuonStart).Mag());
			fmuon_lambda_angle->Fill((LambdaMom.Unit())*(MuonMom.Unit()));
			fproton_pion_angle->Fill((ProtonMom.Unit())*(PionMom.Unit()));
			flambda_KE->Fill((lLambdaMom.E())-(LambdaMass/1000));
			
			
			
			// }}}}}}}}}}
			
			
			// for MCCheater
			
			// grab the sim::ParticleList
			bt->SetEveIdCalculator(new sim::EmEveIdCalculator);
			sim::ParticleList plist = bt->ParticleList();

			// print the list of particles first
			//mf::LogInfo("CheckBackTracking") << plist;
			
			//get the cheated hits
			// grab the hits that have been reconstructed
			art::Handle< std::vector<recob::Hit> > hitcol;
			evt.getByLabel(fHitModuleLabel, hitcol);
			
			// make a vector of them - we aren't writing anything out to a file
			// so no need for a art::PtrVector here
			std::vector< art::Ptr<recob::Hit> > hits;
			art::fill_ptr_vector(hits, hitcol);
			
			art::ServiceHandle<geo::Geometry> geo;
			
			
			// Now get the reco hits
			std::vector< art::Ptr<recob::Hit> > hitlist;
			std::vector< art::Ptr<recob::Hit> > reco_hits;
			
			art::Handle< std::vector<recob::Cluster> > clusterListHandle;
			evt.getByLabel(fLineMergerModuleLabel,clusterListHandle);
			
			art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fLineMergerModuleLabel);

			art::PtrVector<recob::Cluster> cluIn;
			for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii){
			  art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
			  //cluIn.push_back(cluster);
			  if (cluster->View() == geo::kW) continue;
			  cluIn.push_back(cluster);
			  // ptrvector to hits in one cluster
			  hitlist = fmh.at(ii);
			  
			  for(std::vector< art::Ptr<recob::Hit> >::const_iterator theHit = hitlist.begin(); 
			      theHit != hitlist.end();  theHit++) 
			    reco_hits.push_back(*theHit);
			}

			
			double muon_hits_cheated = 0;
			double proton_hits_cheated = 0;
			double pion_hits_cheated = 0;
			
			double muon_hits_cheated_U = 0;
			double proton_hits_cheated_U = 0;
			double pion_hits_cheated_U = 0;
			
			double muon_charge_cheated_U = 0;
			double proton_charge_cheated_U = 0;
			double pion_charge_cheated_U = 0;
			
			double muon_hits_cheated_V = 0;
			double proton_hits_cheated_V = 0;
			double pion_hits_cheated_V = 0;
			
			double muon_charge_cheated_V = 0;
			double proton_charge_cheated_V = 0;
			double pion_charge_cheated_V = 0;
			
			double muon_hits_reco = 0;
			double proton_hits_reco = 0;
			double pion_hits_reco = 0;
			double muon_hits_reco_track = 0;
			double proton_hits_reco_track = 0;
			double pion_hits_reco_track = 0;
			
			
			std::cout << "NUMBER OF CHEATED HITS = " << hits.size() << std::endl;
			
			// loop over the hits and figure out which particle contributed to each one
			std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
			
			unsigned int p(0),w(0), t(0),channel(0), cs(0);
			while( itr != hits.end() ){
			  
			  channel=(*itr)->Wire()->RawDigit()->Channel();
			  geo->ChannelToWire(channel,cs,t,p,w);
			  
			  //std::cout << "____________ HERE IS THE PLANE _____________" << p << std::endl;
			  
			  std::vector<cheat::TrackIDE> eveides   = bt->HitToEveID(*itr);

			  //        std::cout << "EVE ID VECTOR SIZE = " << eveides.size() << std::endl;
			  std::vector<cheat::TrackIDE>::iterator idesitr = eveides.begin();
			  while( idesitr != eveides.end() ){
			    // 	 mf::LogInfo("CheckBackTracking") << "eve id: " << (*idesitr).trackID 
			    // 					  << " contributed " << (*idesitr).energyFrac 
			    // 					  << " to the current hit";
			    
			    //std::cout << "____________ HERE IS THE PLANE _____________" << p << std::endl;
			    
			    // 	 if((*idesitr).trackID == MuonID)
			    // 	   muon_hits_cheated += (*idesitr).energyFrac;
			    
			    // 	 if((*idesitr).trackID == ProtonID)
			    // 	   proton_hits_cheated += (*idesitr).energyFrac;
			    
			    // 	 if((*idesitr).trackID == PionID)
			    // 	   pion_hits_cheated += (*idesitr).energyFrac;
			    
			    
			    if((*idesitr).trackID == MuonID){
			      if(p==0){
				muon_hits_cheated_U += (*idesitr).energyFrac;
				muon_charge_cheated_U += (((*itr)->Charge())*((*idesitr).energyFrac));
				//std::cout << "MUON CHEATED CHARGE IN U = " << muon_charge_cheated_U << std::endl;
			      }
			      if(p==1){
				muon_hits_cheated_V += (*idesitr).energyFrac;
				muon_charge_cheated_V += (((*itr)->Charge())*((*idesitr).energyFrac));
				//std::cout << "MUON CHEATED CHARGE IN V = " << muon_charge_cheated_V << std::endl;
			      }
			    }
			    
			    if((*idesitr).trackID == ProtonID){
			      if(p==0){
				proton_hits_cheated_U += (*idesitr).energyFrac;
				proton_charge_cheated_U += (((*itr)->Charge())*((*idesitr).energyFrac));
			      }
			      if(p==1){
				proton_hits_cheated_V += (*idesitr).energyFrac;
				proton_charge_cheated_V += (((*itr)->Charge())*((*idesitr).energyFrac));
			      }
			    }
			    
			    if((*idesitr).trackID == PionID){
			      if(p==0){
				pion_hits_cheated_U += (*idesitr).energyFrac;
				pion_charge_cheated_U += (((*itr)->Charge())*((*idesitr).energyFrac));
			      }
			      if(p==1){
				pion_hits_cheated_V += (*idesitr).energyFrac;
				pion_charge_cheated_V += (((*itr)->Charge())*((*idesitr).energyFrac));
			      }
			    }
			    
			    idesitr++;
			    
			  }
			  
			  itr++;
			  
			}
			
			muon_hits_cheated = muon_hits_cheated_U + muon_hits_cheated_V;
			proton_hits_cheated = proton_hits_cheated_U + proton_hits_cheated_V;
			pion_hits_cheated = pion_hits_cheated_U + pion_hits_cheated_V;
			
			std::cout << "NUMBER OF RECO HITS = " << reco_hits.size() << std::endl;
			
			// loop over the all the hits from merged lines and figure out which particle contributed to each hit
			std::vector< art::Ptr<recob::Hit> >::iterator reco_itr = reco_hits.begin();
			
			
			while( reco_itr != reco_hits.end() ){
			  
			  //std::cout << "I AM GOING TO CALCULATE THE RECO_EVE_IDS" << std::endl;
			  std::vector<cheat::TrackIDE> reco_eveides   = bt->HitToEveID(*reco_itr);
			  //        std::cout << "EVE ID VECTOR SIZE RECO = " << reco_eveides.size() << std::endl;
			  std::vector<cheat::TrackIDE>::iterator reco_idesitr = reco_eveides.begin();
			  while( reco_idesitr != reco_eveides.end() ){
			    //std::cout << "I AM GOING TO PRINT THE RECO_EVE_IDS" << std::endl;
			    // 	 mf::LogInfo("CheckBackTracking") << "eve id: " << (*reco_idesitr).trackID 
			    // 					  << " contributed " << (*reco_idesitr).energyFrac 
			    // 					  << " to the current hit";
			    
			    if((*reco_idesitr).trackID == MuonID)
			      muon_hits_reco += (*reco_idesitr).energyFrac ;
			    
			    if((*reco_idesitr).trackID == ProtonID)
			      proton_hits_reco += (*reco_idesitr).energyFrac ;
			    
			    if((*reco_idesitr).trackID == PionID)
			      pion_hits_reco += (*reco_idesitr).energyFrac ;
			    
			    reco_idesitr++;
			    
			  }
			  
			  reco_itr++;
			  
			}
			
			// Now see the track_reco  hits
			std::vector< art::Ptr<recob::Hit> > track_reco_hits;
			art::PtrVector<recob::Cluster> Icluster;
			art::PtrVector<recob::Cluster> Ccluster;
			
			for(size_t i = 0; i < trkIn.size(); ++i){
			  std::vector< art::Ptr<recob::Hit> > track_hitlist = fmht.at(i);
			  
			  for(size_t h = 0; h < track_hitlist.size(); ++h)
			    track_reco_hits.push_back(track_hitlist[h]);
			}
			
			std::cout << "NUMBER OF track_RECO HITS = " << track_reco_hits.size() << std::endl;
			
			// loop over the hits and figure out which particle contributed to each one
			std::vector< art::Ptr<recob::Hit> >::iterator track_reco_itr = track_reco_hits.begin();
			
			
			while( track_reco_itr != track_reco_hits.end() ){
			  //std::cout << "I AM GOING TO CALCULATE THE RECO_EVE_IDS" << std::endl;
			  std::vector<cheat::TrackIDE> track_reco_eveides   = bt->HitToEveID(*track_reco_itr);
			  //        std::cout << "EVE ID VECTOR SIZE RECO = " << reco_eveides.size() << std::endl;
			  std::vector<cheat::TrackIDE>::iterator track_reco_idesitr = track_reco_eveides.begin();
			  while( track_reco_idesitr != track_reco_eveides.end() ){
			    //std::cout << "I AM GOING TO PRINT THE RECO_EVE_IDS" << std::endl;
			    // 	 mf::LogInfo("CheckBackTracking") << "eve id: " << (*reco_idesitr).trackID 
			    // 					  << " contributed " << (*reco_idesitr).energyFrac 
			    // 					  << " to the current hit";
			    
			    if((*track_reco_idesitr).trackID == MuonID)
			      muon_hits_reco_track += (*track_reco_idesitr).energyFrac ;
			    
			    if((*track_reco_idesitr).trackID == ProtonID)
			      proton_hits_reco_track += (*track_reco_idesitr).energyFrac ;
			    
			    if((*track_reco_idesitr).trackID == PionID)
			      pion_hits_reco_track += (*track_reco_idesitr).energyFrac ;
			    
			    track_reco_idesitr++;
			    
			  }
			  
			  track_reco_itr++;
			  
			}
			
			
			//plot hists
			
			int uclus=0; 
			int vclus =0;
			art::PtrVector<recob::Cluster> ucluIn;
			art::PtrVector<recob::Cluster> vcluIn;
			
			for(unsigned int i=0; i<cluIn.size();++i){
			  
			  switch(cluIn[i]->View()){
			  case geo::kU :
			    uclus ++;
			    ucluIn.push_back(cluIn[i]);
			    break;
			  case geo::kV :
			    vclus ++;
			    vcluIn.push_back(cluIn[i]);
			    break;
			  default :
			    break;
			  }
			}
			
			art::FindManyP<recob::Hit> fmhuc(ucluIn, evt, fLineMergerModuleLabel);
			art::FindManyP<recob::Hit> fmhvc(vcluIn, evt, fLineMergerModuleLabel);

			std::cout << "number of clusters in U view = " << uclus << std::endl;
			std::cout << "number of clusters in V view = " << vclus << std::endl;
			
			fU_no_clus->Fill(uclus);
			fV_no_clus->Fill(vclus);
			
			
			std::vector<int> skip_mu_U;
			std::vector<int> skip_p_U;
			std::vector<int> skip_pi_U;
			std::vector<int> skip_mu_V;
			std::vector<int> skip_p_V;
			std::vector<int> skip_pi_V;
			
			double selected_mu_hits_U = 0.;
			double selected_p_hits_U = 0.;
			double selected_pi_hits_U = 0.;
			double selected_mu_hits_V = 0.;
			double selected_p_hits_V = 0.;
			double selected_pi_hits_V = 0.;
			
			double selected_mu_charge_U = 0.;
			double selected_p_charge_U = 0.;
			double selected_pi_charge_U = 0.;
			double selected_mu_charge_V = 0.;
			double selected_p_charge_V = 0.;
			double selected_pi_charge_V = 0.;
			
			for(unsigned int  i = 0; i < ucluIn.size();  i++){
			  double startwire = ucluIn[i]->StartPos()[0];
			  double starttime = ucluIn[i]->StartPos()[1];
			  double endwire = ucluIn[i]->EndPos()[0];
			  double endtime = ucluIn[i]->EndPos()[1];
			  double uclu_length = sqrt((pow(startwire-endwire,2)*13.5)+pow(starttime-endtime,2));
			  //std::cout << "ucluster_length = " << uclu_length << std::endl;
			  fU_clu_length->Fill(uclu_length);
			  fUwire_span_length->Fill(abs(endwire-startwire)*0.4);
			  
			  std::vector< art::Ptr<recob::Hit> > uclu_hits;
			  std::vector< art::Ptr<recob::Hit> > reco_hits_U;
			  double total_charge_U = 0.;
			  int total_number_of_hits_U = 0.;
			  double mu_hits_U = 0.;
			  double mu_charge_U = 0.;
			  double p_hits_U = 0.;
			  double p_charge_U = 0.;
			  double pi_hits_U = 0.;
			  double pi_charge_U = 0.;
			  
			  //calculate purity
			  uclu_hits = fmhuc.at(i);
			  for(std::vector< art::Ptr<recob::Hit> >::const_iterator theHit = uclu_hits.begin(); theHit != uclu_hits.end();  theHit++){
			    reco_hits_U.push_back(*theHit);
			  }
			  
			  for(unsigned int  j= 0; j< reco_hits_U.size(); ++j){
			    //std::vector< art::Ptr<recob::Hit> >::iterator reco_itr_U = reco_hits_U.begin();
			    //while( reco_itr_U != reco_hits_U.end() ){
			    
			    total_charge_U += reco_hits_U[j]->Charge();
			    total_number_of_hits_U = reco_hits_U.size();
			    
			    std::vector<cheat::TrackIDE> reco_eveides_U   =bt->HitToEveID(reco_hits_U[j]);
			    std::vector<cheat::TrackIDE>::iterator reco_idesitr_U = reco_eveides_U.begin();
			    
			    //std::cout<< "eve id size = " << reco_eveides_U.size() << std::endl;
			    
			    while(reco_idesitr_U != reco_eveides_U.end()){
			      
			      //  	 mf::LogInfo("CheckBackTracking") << "eve id: " << (*reco_idesitr_U).trackID 
			      //  					  << " contributed " << (*reco_idesitr_U).energyFrac; 
			      
			      
			      if((*reco_idesitr_U).trackID == MuonID){
				mu_hits_U += (*reco_idesitr_U).energyFrac ;
				mu_charge_U += ((reco_hits_U[j]->Charge()) * ((*reco_idesitr_U).energyFrac)) ;
			      }
			      
			      if((*reco_idesitr_U).trackID == ProtonID){
				p_hits_U += (*reco_idesitr_U).energyFrac ;
				p_charge_U += ((reco_hits_U[j]->Charge()) * ((*reco_idesitr_U).energyFrac)) ;
			      }
			      
			      if((*reco_idesitr_U).trackID == PionID){
				pi_hits_U += (*reco_idesitr_U).energyFrac ;
				pi_charge_U += ((reco_hits_U[j]->Charge()) * ((*reco_idesitr_U).energyFrac)) ;
			      }
			      
			      reco_idesitr_U++;
			      
			    }//while
			    
			  }//loop over hits
			  
			  std::cout << "***********************" << std::endl;
			  std::cout << "total charge of this Icluster = " << total_charge_U << std::endl;
			  std::cout << "total, mu , p and pi hits are "<< total_number_of_hits_U 
				    << " " << mu_hits_U << " "<< p_hits_U<< " " << pi_hits_U << std::endl;
			  
			  if((mu_hits_U/total_number_of_hits_U)>=0.9){
			    fpurity_muon->Fill((lMuonMom.E()-(MuonMass/1000)), (mu_charge_U/total_charge_U));
			    fcompleteness_muon->Fill((lMuonMom.E()-(MuonMass/1000)), (mu_charge_U/muon_charge_cheated_U));
			    std::cout<< "ABOVE: I AM A MUON IN INDUCTION" << std::endl;
			    if(skip_mu_U.size()==0){
			      std::cout<< "ABOVE: I AM A SELECTED MUON IN INDUCTION" << std::endl;
			      selected_mu_hits_U = mu_hits_U;
			      selected_mu_charge_U = total_charge_U;
			      skip_mu_U.push_back(1);
			    }
			  }
			  else if((p_hits_U/total_number_of_hits_U)>=0.9){
			    fpurity_proton->Fill((lProtonMom.E()-(ProtonMass/1000)), (p_charge_U/total_charge_U));
			    fcompleteness_proton->Fill((lProtonMom.E()-(ProtonMass/1000)), (p_charge_U/proton_charge_cheated_U));
			    std::cout<< "ABOVE: I AM A PROTON IN INDUCTION" << std::endl;
			    if(skip_p_U.size()==0){
			      std::cout<< "ABOVE: I AM A SELECTED PROTON IN INDUCTION" << std::endl;
			      selected_p_hits_U = p_hits_U;
			      selected_p_charge_U = total_charge_U;
			      skip_p_U.push_back(1);
			    }
			  }
			  else if((pi_hits_U/total_number_of_hits_U)>=0.9){
			    fpurity_pion->Fill((lPionMom.E()-(PionMass/1000)), (pi_charge_U/total_charge_U));
			    fcompleteness_pion->Fill((lPionMom.E()-(PionMass/1000)), (pi_charge_U/pion_charge_cheated_U));
			    std::cout<< "ABOVE: I AM A PION IN INDUCTION" << std::endl;
			    if(skip_pi_U.size()==0){
			      std::cout<< "ABOVE: I AM A SELECTED PION IN INDUCTION" << std::endl;
			      selected_pi_hits_U = pi_hits_U;
			      selected_pi_charge_U = total_charge_U;
			      skip_pi_U.push_back(1);
			    }
			  }
			  
			}//loop over u clusters
			
			

			
			
			for(unsigned int  i = 0; i < vcluIn.size();  i++){
			  double startwire = vcluIn[i]->StartPos()[0];
			  double starttime = vcluIn[i]->StartPos()[1];
			  double endwire = vcluIn[i]->EndPos()[0];
			  double endtime = vcluIn[i]->EndPos()[1];
			  double vclu_length = sqrt((pow(startwire-endwire,2)*13.5)+pow(starttime-endtime,2));
			  //std::cout << "vcluster_length = " << vclu_length << std::endl;
			  fV_clu_length->Fill(vclu_length);
			  fVwire_span_length->Fill(abs(endwire-startwire)*0.4);
			  
			  std::vector< art::Ptr<recob::Hit> > vclu_hits;
			  std::vector< art::Ptr<recob::Hit> > reco_hits_V;
			  double total_charge_V = 0.;
			  int total_number_of_hits_V = 0.;
			  double mu_hits_V = 0.;
			  double mu_charge_V = 0.;
			  double p_hits_V = 0.;
			  double p_charge_V = 0.;
			  double pi_hits_V = 0.;
			  double pi_charge_V = 0.;
			  
			  //calculate purity
			  vclu_hits = fmhvc.at(i);
			  for(std::vector< art::Ptr<recob::Hit> >::const_iterator theHit = vclu_hits.begin(); theHit != vclu_hits.end();  theHit++){
			    reco_hits_V.push_back(*theHit);
			  }
			  
			  for(unsigned int  j= 0; j< reco_hits_V.size(); ++j){
			    //std::vector< art::Ptr<recob::Hit> >::iterator reco_itr_U = reco_hits_U.begin();
			    //while( reco_itr_U != reco_hits_U.end() ){
			    
			    total_charge_V += reco_hits_V[j]->Charge();
			    total_number_of_hits_V = reco_hits_V.size();
			    
			    std::vector<cheat::TrackIDE> reco_eveides_V   = bt->HitToEveID(reco_hits_V[j]);
			    std::vector<cheat::TrackIDE>::iterator reco_idesitr_V = reco_eveides_V.begin();
			    
			    //std::cout<< "eve id size = " << reco_eveides_V.size() << std::endl;
			    
			    while(reco_idesitr_V != reco_eveides_V.end()){
			      
			      //  	 mf::LogInfo("CheckBackTracking") << "eve id: " << (*reco_idesitr_V).trackID 
			      //  					  << " contributed " << (*reco_idesitr_V).energyFrac; 
			      
			      
			      if((*reco_idesitr_V).trackID == MuonID){
				mu_hits_V += (*reco_idesitr_V).energyFrac ;
				mu_charge_V += ((reco_hits_V[j]->Charge()) * ((*reco_idesitr_V).energyFrac)) ;
			      }
			      
			      if((*reco_idesitr_V).trackID == ProtonID){
				p_hits_V += (*reco_idesitr_V).energyFrac ;
				p_charge_V += ((reco_hits_V[j]->Charge()) * ((*reco_idesitr_V).energyFrac)) ;
			      }
			      
			      if((*reco_idesitr_V).trackID == PionID){
				pi_hits_V += (*reco_idesitr_V).energyFrac ;
				pi_charge_V += ((reco_hits_V[j]->Charge()) * ((*reco_idesitr_V).energyFrac)) ;
			      }
			      
			      reco_idesitr_V++;
			      
			    }//while
			    
			  }//loop over hits
			  
			  std::cout << "***********************" << std::endl;
			  std::cout << "total charge of this Ccluster = " << total_charge_V << std::endl;
			  std::cout << "total, mu , p and pi hits are "<< total_number_of_hits_V << " " << mu_hits_V << " "<< p_hits_V<< " " << pi_hits_V << std::endl;
			  
			  if((mu_hits_V/total_number_of_hits_V)>=0.9){
			    fpurity_muon->Fill((lMuonMom.E()-(MuonMass/1000)), (mu_charge_V/total_charge_V));
			    fcompleteness_muon->Fill((lMuonMom.E()-(MuonMass/1000)), (mu_charge_V/muon_charge_cheated_V)); 
			    std::cout<< "ABOVE: I AM A MUON IN COLLECTION" << std::endl;
			    if(skip_mu_V.size()==0){
			      std::cout<< "ABOVE: I AM A SELECTED MUON IN COLLECTION" << std::endl;
			      selected_mu_hits_V = mu_hits_V;
			      selected_mu_charge_V = total_charge_V;
			      skip_mu_V.push_back(1);
			    }
			  }
			  else if((p_hits_V/total_number_of_hits_V)>=0.9){
			    fpurity_proton->Fill((lProtonMom.E()-(ProtonMass/1000)), (p_charge_V/total_charge_V));
			    fcompleteness_proton->Fill((lProtonMom.E()-(ProtonMass/1000)), (p_charge_V/proton_charge_cheated_V));
			    std::cout<< "ABOVE: I AM A PROTON IN COLLECTION" << std::endl;
			    if(skip_p_V.size()==0){
			      std::cout<< "ABOVE: I AM A SELECTED PROTON IN COLLECTION" << std::endl;
			      selected_p_hits_V = p_hits_V;
			      selected_p_charge_V = total_charge_V;
			      skip_p_V.push_back(1);
			    }
			  }
			  else if((pi_hits_V/total_number_of_hits_V)>=0.9){
			    fpurity_pion->Fill((lPionMom.E()-(PionMass/1000)), (pi_charge_V/total_charge_V));
			    fcompleteness_pion->Fill((lPionMom.E()-(PionMass/1000)), (pi_charge_V/pion_charge_cheated_V));
			    std::cout<< "ABOVE: I AM A PION IN COLLECTION" << std::endl;
			    if(skip_pi_V.size()==0){
			      std::cout<< "ABOVE: I AM A SELECTED PION IN COLLECTION" << std::endl;
			      selected_pi_hits_V = pi_hits_V;
			      selected_pi_charge_V = total_charge_V;
			      skip_pi_V.push_back(1);
			    }
			  }
			}//loop over v clusters
			
			
			std::cout << "___________________________" << std::endl;
			std::cout << "mu hits in U = " << selected_mu_hits_U << " mu hits in V = " << selected_mu_hits_V << std::endl;
			std::cout << "p hits in U = " << selected_p_hits_U << " p hits in V = " << selected_p_hits_V << std::endl;
			std::cout << "pi hits in U = " << selected_pi_hits_U << " pi hits in V = " << selected_pi_hits_V << std::endl;
			std::cout << "___________________________" << std::endl;
			
			std::cout << "mu charge in U = " << selected_mu_charge_U << " mu charge in V = " << selected_mu_charge_V << std::endl;
			std::cout << "p charge in U = " << selected_p_charge_U << " p charge in V = " << selected_p_charge_V << std::endl;
			std::cout << "pi charge in U = " << selected_pi_charge_U << " pi charge in V = " << selected_pi_charge_V << std::endl;
			std::cout << "___________________________" << std::endl;
			
			
			fmuon_cluster_hits_UV->Fill(selected_mu_hits_U,selected_mu_hits_V);
			fproton_cluster_hits_UV->Fill(selected_p_hits_U,selected_p_hits_V);
			fpion_cluster_hits_UV->Fill(selected_pi_hits_U,selected_pi_hits_V);
			
			fmuon_cluster_charge_UV->Fill((selected_mu_charge_U/1000),(selected_mu_charge_V/1000));
			fproton_cluster_charge_UV->Fill((selected_p_charge_U/1000),(selected_p_charge_V/1000));
			fpion_cluster_charge_UV->Fill((selected_pi_charge_U/1000),(selected_pi_charge_V/1000));
			//     TProfile *prof = fmuon_cluster_charge_UV->ProfileX();
			//     prof->Fit("pol1");
			
			if(selected_mu_hits_U==0 || selected_mu_hits_V==0){
			  fmuon_reco_percentage->Fill(0);
			}
			else
			  fmuon_reco_percentage->Fill(1);
			
			if(selected_p_hits_U==0 || selected_p_hits_V==0){
			  fproton_reco_percentage->Fill(0);
			}
			else
			  fproton_reco_percentage->Fill(1);
			
			if(selected_pi_hits_U==0 || selected_pi_hits_V==0){
			  fpion_reco_percentage->Fill(0);
			}
			else
			  fpion_reco_percentage->Fill(1);
			
			
			
			std::cout << "muon_hits_cheated = " << muon_hits_cheated << std::endl;
			std::cout << "muon_hits_reco = " << muon_hits_reco << std::endl;
			std::cout << "muon_hits_reco_track = " << muon_hits_reco_track << std::endl;
			
			std::cout << "muon_hits_cheated_U = " << muon_hits_cheated_U << std::endl;
			std::cout << "muon_hits_cheated_V = " << muon_hits_cheated_V << std::endl;
			
			std::cout << "proton_hits_cheated = " << proton_hits_cheated << std::endl;
			std::cout << "proton_hits_reco = " << proton_hits_reco << std::endl;
			std::cout << "proton_hits_reco_track = " << proton_hits_reco_track << std::endl;
			
			std::cout << "pion_hits_cheated = " << pion_hits_cheated << std::endl;
			std::cout << "pion_hits_reco = " << pion_hits_reco << std::endl;
			std::cout << "pion_hits_reco_track = " << pion_hits_reco << std::endl;
			
			std::cout << "number of 3d tracks = " << trkIn.size() << std::endl;
			
			
			fMuon_Reco_Eff->Fill((lMuonMom.E()-(MuonMass/1000)),(muon_hits_reco/muon_hits_cheated));
			fProtonn_Reco_Eff->Fill((lProtonMom.E()-(ProtonMass/1000)),(proton_hits_reco/proton_hits_cheated));
			fPion_Reco_Eff->Fill((lPionMom.E()-(PionMass/1000)),(pion_hits_reco/pion_hits_cheated));
			
			fMuon_Reco_Eff_vs_hits->Fill(muon_hits_cheated,(muon_hits_reco/muon_hits_cheated));
			fProtonn_Reco_Eff_vs_hits->Fill(proton_hits_cheated,(proton_hits_reco/proton_hits_cheated));
			fPion_Reco_Eff_vs_hits->Fill(pion_hits_cheated,(pion_hits_reco/pion_hits_cheated));
			
			fMuon_Reco_Eff_3d->Fill((lMuonMom.E()-(MuonMass/1000)),(muon_hits_reco_track/muon_hits_cheated));
			fProtonn_Reco_Eff_3d->Fill((lProtonMom.E()-(ProtonMass/1000)),(proton_hits_reco_track/proton_hits_cheated));
			fPion_Reco_Eff_3d->Fill((lPionMom.E()-(PionMass/1000)),(pion_hits_reco_track/pion_hits_cheated));
			
			fNumber_tracks->Fill(trkIn.size());
			
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

		
    
  } // analyze
} // namespace

////////////////////////////////////////////////////////////////////////
