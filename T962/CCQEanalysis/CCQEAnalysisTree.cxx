////////////////////////////////////////////////////////////////////////
//
// Create a TTree for analysis
//
// \author tjyang@fnal.gov
// \author joshua.spitz@yale.edu
// \author kinga.partyka@yale.edu
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


#include "T962/CCQEanalysis/CCQEAnalysisTree.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "SummaryData/summary.h"

 
//-------------------------------------------------
t962::CCQEAnalysisTree::CCQEAnalysisTree(fhicl::ParameterSet const& pset) : 
  
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fKingaModuleLabel         (pset.get< std::string >("KingaModuleLabel")        ),
  fLineMergerModuleLabel       (pset.get< std::string >("LineMergerModuleLabel")      ),
  fDbscanModuleLabel       (pset.get< std::string >("DbscanModuleLabel")  
  ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fMINOSModuleLabel         (pset.get< std::string >("MINOSModuleLabel")        ),
  fTrackMatchModuleLabel    (pset.get< std::string >("TrackMatchModuleLabel")   ),
  fScanModuleLabel          (pset.get< std::string >("ScanModuleLabel")         ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fvertextrackWindow        (pset.get< double >("vertextrackWindow")            ),
  fvertexclusterWindow      (pset.get< double >("vertexclusterWindow")          ),
  fboundaryWindow           (pset.get< double >("boundaryWindow")               ),
  no_kingaclusters(400),
  no_linemergerclusters(400),
  ntracks_reco(400),
  no_primaries(400),
  genie_no_primaries(400),
  no_hits(20000)
 
  
{
}

//-------------------------------------------------
t962::CCQEAnalysisTree::~CCQEAnalysisTree()
{
delete twodvtx_w_reco;
delete twodvtx_t_reco;
delete twodvtx_w_truth;
delete twodvtx_t_truth;
delete Start_pt_w_kingaCl;
delete Start_pt_t_kingaCl;
delete kingaclusters_planeNo;
delete linemergerclusters_planeNo;
delete Start_pt_w_linemergerCl;
delete Start_pt_t_linemergerCl;
delete two_trackstart_dcosx_reco;
delete two_trackstart_dcosy_reco;
delete two_trackstart_dcosz_reco;
delete two_trackexit_dcosx_reco;
delete two_trackexit_dcosy_reco;
delete two_trackexit_dcosz_reco;
delete all_trackstart_x_reco;
delete all_trackstart_y_reco;
delete all_trackstart_z_reco;
delete primaries_pdg;
delete Eng;
delete Px;
 delete Py;
 delete Pz;
 delete StartPointx;
 delete StartPointy;
 delete StartPointz;
 delete EndPointx;
 delete EndPointy;
 delete EndPointz;
 delete NumberDaughters;
 delete genie_primaries_pdg;
delete genie_Eng;
delete genie_Px;
 delete genie_Py;
 delete genie_Pz;
 delete genie_P;
 delete genie_status_code;
 delete genie_mass;
 delete genie_trackID;
 delete genie_ND;
 delete genie_mother;
 delete hit_plane;
 delete hit_wire;
 delete hit_channel;
 delete hit_peakT;
 delete hit_charge;
}

void t962::CCQEAnalysisTree::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  
  all_trackstart_x_reco=new double[ntracks_reco];
  all_trackstart_y_reco=new double[ntracks_reco];
  all_trackstart_z_reco=new double[ntracks_reco];
  twodvtx_w_reco= new double[2];
  twodvtx_t_reco= new double[2];
  twodvtx_w_truth= new double[2];
  twodvtx_t_truth= new double[2];
  kingaclusters_planeNo=new int[no_kingaclusters];
  linemergerclusters_planeNo=new int[no_linemergerclusters];
  Start_pt_w_kingaCl=new double[no_kingaclusters];
  Start_pt_t_kingaCl=new double[no_kingaclusters];
  Start_pt_w_linemergerCl=new double[no_linemergerclusters];
  Start_pt_t_linemergerCl=new double[no_linemergerclusters];
  two_trackstart_dcosx_reco= new double[2];
  two_trackstart_dcosy_reco= new double[2];
  two_trackstart_dcosz_reco= new double[2];
  two_trackexit_dcosx_reco= new double[2];
  two_trackexit_dcosy_reco= new double[2];
  two_trackexit_dcosz_reco= new double[2];
  
   primaries_pdg= new int[no_primaries];
  Eng= new double[no_primaries];
  Px= new double[no_primaries];
  Py= new double[no_primaries];
  Pz= new double[no_primaries];
  StartPointx= new double[no_primaries];
  StartPointy= new double[no_primaries];
  StartPointz= new double[no_primaries];
  EndPointx= new double[no_primaries];
  EndPointy= new double[no_primaries];
  EndPointz= new double[no_primaries];
  NumberDaughters= new int[no_primaries];
 
 
genie_primaries_pdg= new double[genie_no_primaries];
 genie_Eng= new double[genie_no_primaries];
 genie_Px= new double[genie_no_primaries];
  genie_Py= new double[genie_no_primaries];
 genie_Pz= new double[genie_no_primaries];
 genie_P= new double[genie_no_primaries];
genie_status_code= new int[genie_no_primaries];
genie_mass= new double[genie_no_primaries];
genie_trackID= new int[genie_no_primaries];
genie_ND= new int[genie_no_primaries];
 genie_mother= new int[genie_no_primaries];
 
 hit_plane= new int[no_hits];
   hit_wire= new int[no_hits];
    hit_channel= new int[no_hits];
   hit_peakT= new double[no_hits];
   hit_charge= new double[no_hits];
   
    
  fTree->Branch("no_hits",&no_hits,"no_hits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[no_hits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[no_hits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[no_hits]/I");
   fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/D");
   fTree->Branch("hit_charge",hit_charge,"hit_charge[no_hits]/D");
  
  
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("vtxx_reco",&vtxx_reco,"vtxx_reco/D");
  fTree->Branch("vtxy_reco",&vtxy_reco,"vtxy_reco/D");
  fTree->Branch("vtxz_reco",&vtxz_reco,"vtxz_reco/D");
  
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("nvertextracks_reco",&nvertextracks_reco,"nvertextracks_reco/I");
 
  fTree->Branch("ntrackendonboundary_reco",&ntrackendonboundary_reco,"ntrackendonboundary_reco/I");
  fTree->Branch("trackstart_dcosx_reco",&trackstart_dcosx_reco, "trackstart_dcosx_reco/D");
  fTree->Branch("trackstart_dcosy_reco",&trackstart_dcosy_reco, "trackstart_dcosy_reco/D");
  fTree->Branch("trackstart_dcosz_reco",&trackstart_dcosz_reco, "trackstart_dcosz_reco/D");
  fTree->Branch("trackexit_dcosx_reco",&trackexit_dcosx_reco, "trackexit_dcosx_reco/D");
  fTree->Branch("trackexit_dcosy_reco",&trackexit_dcosy_reco, "trackexit_dcosy_reco/D");
  fTree->Branch("trackexit_dcosz_reco",&trackexit_dcosz_reco, "trackexit_dcosz_reco/D");
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
  
  
  fTree->Branch("vtxx_scan", &vtxx_scan, "vtxx_scan/D");
  fTree->Branch("vtxy_scan", &vtxy_scan, "vtxy_scan/D");
  fTree->Branch("vtxz_scan", &vtxz_scan, "vtxz_scan/D");
  fTree->Branch("ntracks_scan", &ntracks_scan , "ntracks_scan/I");
  fTree->Branch("nshowers_scan", &nshowers_scan , "nshowers_scan/I");
  fTree->Branch("neutrino_scan", &neutrino_scan , "neutrino_scan/I");
  fTree->Branch("maybeneutrino_scan", &maybeneutrino_scan , "maybeneutrino_scan/I");    
  fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  fTree->Branch("enu_truth",&enu_truth,"enu_truth/D");
  fTree->Branch("Q2_truth",&Q2_truth,"Q2_truth/D");
  fTree->Branch("W_truth",&W_truth,"W_truth/D");
  fTree->Branch("nuvtxx_truth",&nuvtxx_truth,"nuvtxx_truth/D");
  fTree->Branch("nuvtxy_truth",&nuvtxy_truth,"nuvtxy_truth/D");
  fTree->Branch("nuvtxz_truth",&nuvtxz_truth,"nuvtxz_truth/D");
  fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
  fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
  fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
  fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");
  
  //kinga:
  fTree->Branch("twodvtx_w_reco", twodvtx_w_reco, "twodvtx_w_reco[2]/D");
  fTree->Branch("twodvtx_t_reco", twodvtx_t_reco, "twodvtx_t_reco[2]/D");
  fTree->Branch("twodvtx_w_truth", twodvtx_w_truth, "twodvtx_w_truth[2]/D");
  fTree->Branch("twodvtx_t_truth", twodvtx_t_truth, "twodvtx_t_truth[2]/D");
  fTree->Branch("nkingaclustersu_reco",&nkingaclustersu_reco,"nkingaclustersu_reco/I");
  fTree->Branch("nkingaclustersv_reco",&nkingaclustersv_reco,"nkingaclustersv_reco/I");
  fTree->Branch("nvertexkingaclustersu_reco",&nvertexkingaclustersu_reco,"nvertexkingaclustersu_reco/I");
  fTree->Branch("nvertexkingaclustersv_reco",&nvertexkingaclustersv_reco,"nvertexkingaclustersv_reco/I");
  
  
  fTree->Branch("nlinemergerclustersu_reco",&nlinemergerclustersu_reco,"nlinemergerclustersu_reco/I");
  fTree->Branch("nlinemergerclustersv_reco",&nlinemergerclustersv_reco,"nlinemergerclustersv_reco/I");
  fTree->Branch("nvertexlinemergerclustersu_reco",&nvertexlinemergerclustersu_reco,"nvertexlinemergerclustersu_reco/I");
  fTree->Branch("nvertexlinemergerclustersv_reco",&nvertexlinemergerclustersv_reco,"nvertexlinemergerclustersv_reco/I");
  
  
  fTree->Branch("ndbscanclustersu_reco",&ndbscanclustersu_reco,"ndbscanclustersu_reco/I");
  fTree->Branch("ndbscanclustersv_reco",&ndbscanclustersv_reco,"ndbscanclustersv_reco/I");
  fTree->Branch("nvertexdbscanclustersu_reco",&nvertexdbscanclustersu_reco,"nvertexdbscanclustersu_reco/I");
  fTree->Branch("nvertexdbscanclustersv_reco",&nvertexdbscanclustersv_reco,"nvertexdbscanclustersv_reco/I");
  
  
  
  fTree->Branch("no_kingaclusters",&no_kingaclusters,"no_kingaclusters/I");
  fTree->Branch("kingaclusters_planeNo",kingaclusters_planeNo,"kingaclusters_planeNo[no_kingaclusters]/I");
  fTree->Branch("Start_pt_w_kingaCl", Start_pt_w_kingaCl, "Start_pt_w_kingaCl[no_kingaclusters]/D");
 fTree->Branch("Start_pt_t_kingaCl", Start_pt_t_kingaCl, "Start_pt_t_kingaCl[no_kingaclusters]/D");
 
 fTree->Branch("no_linemergerclusters",&no_linemergerclusters,"no_linemergerclusters/I");
 fTree->Branch("linemergerclusters_planeNo",linemergerclusters_planeNo,"linemergerclusters_planeNo[no_linemergerclusters]/I");
  fTree->Branch("Start_pt_w_linemergerCl", Start_pt_w_linemergerCl, "Start_pt_w_linemergerCl[no_linemergerclusters]/D");
  fTree->Branch("Start_pt_t_linemergerCl", Start_pt_t_linemergerCl, "Start_pt_t_linemergerCl[no_linemergerclusters]/D");
  
  fTree->Branch("two_trackstart_dcosx_reco",two_trackstart_dcosx_reco, "two_trackstart_dcosx_reco[2]/D");
  fTree->Branch("two_trackstart_dcosy_reco",two_trackstart_dcosy_reco, "two_trackstart_dcosy_reco[2]/D");
  fTree->Branch("two_trackstart_dcosz_reco",two_trackstart_dcosz_reco, "two_trackstart_dcosz_reco[2]/D");
  fTree->Branch("two_trackexit_dcosx_reco",two_trackexit_dcosx_reco, "two_trackexit_dcosx_reco[2]/D");
  fTree->Branch("two_trackexit_dcosy_reco",two_trackexit_dcosy_reco, "two_trackexit_dcosy_reco[2]/D");
   fTree->Branch("two_trackexit_dcosz_reco",two_trackexit_dcosz_reco, "two_trackexit_dcosz_reco[2]/D");
 
 fTree->Branch("all_trackstart_x_reco", all_trackstart_x_reco, "all_trackstart_x_reco[ntracks_reco]/D");
  fTree->Branch("all_trackstart_y_reco", all_trackstart_y_reco, "all_trackstart_y_reco[ntracks_reco]/D");
  fTree->Branch("all_trackstart_z_reco", all_trackstart_z_reco, "all_trackstart_z_reco[ntracks_reco]/D");
  
    //......................................................
// from geant4:

  fTree->Branch("no_primaries",&no_primaries,"no_primaries/I");
  fTree->Branch("primaries_pdg",primaries_pdg,"primaries_pdg[no_primaries]/I");
  fTree->Branch("Eng",Eng,"Eng[no_primaries]/D");
  fTree->Branch("Px",Px,"Px[no_primaries]/D");
  fTree->Branch("Py",Py,"Py[no_primaries]/D");
  fTree->Branch("Pz",Pz,"Pz[no_primaries]/D");
  fTree->Branch("StartPointx",StartPointx,"StartPointx[no_primaries]/D");
  fTree->Branch("StartPointy",StartPointy,"StartPointy[no_primaries]/D");
  fTree->Branch("StartPointz",StartPointz,"StartPointz[no_primaries]/D");
  fTree->Branch("EndPointx",EndPointx,"EndPointx[no_primaries]/D");
  fTree->Branch("EndPointy",EndPointy,"EndPointy[no_primaries]/D");
  fTree->Branch("EndPointz",EndPointz,"EndPointz[no_primaries]/D");
  fTree->Branch("NumberDaughters",NumberDaughters,"NumberDaughters[no_primaries]/I");
  
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  
  
 //..................................
 // now from genie:
  fTree->Branch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
  fTree->Branch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/D");
  fTree->Branch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/D");
  fTree->Branch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/D");
  fTree->Branch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/D");
  fTree->Branch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/D");
  fTree->Branch("genie_P",genie_P,"genie_P[genie_no_primaries]/D");
  fTree->Branch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/D");
  fTree->Branch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/D");
  fTree->Branch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
  fTree->Branch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
  fTree->Branch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");
  
 
 
}


void t962::CCQEAnalysisTree::beginSubRun(const art::SubRun& sr)
{

art::Handle< sumdata::POTSummary > potListHandle;
sr.getByLabel(fPOTModuleLabel,potListHandle);

if(sr.getByLabel(fPOTModuleLabel,potListHandle))
pot=potListHandle->totpot;
else
pot=0.;

}




void t962::CCQEAnalysisTree::analyze(const art::Event& evt)
{
std::cout<<" IN *** MY *** CCQEANALYSISTREE ***"<<std::endl;
  ResetVars();

  run = evt.run();
  event = evt.id().event();

  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;

  art::Handle< std::vector<raw::RawDigit> > rdListHandle;
  evt.getByLabel(fDigitModuleLabel,rdListHandle);
  art::Handle< std::vector<sim::SimChannel> > scListHandle;
  evt.getByLabel(fDigitModuleLabel,scListHandle);
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  std::vector< art::Ptr<recob::Hit> > hits;
  art::fill_ptr_vector(hits, hitListHandle);
  
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  art::Handle< std::vector<recob::Wire> > wireListHandle;
  evt.getByLabel(fCalDataModuleLabel,wireListHandle);
  art::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrackModuleLabel,trackListHandle);
  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  art::Handle< std::vector<t962::MINOS> > minosListHandle;
  evt.getByLabel(fMINOSModuleLabel,minosListHandle);
  art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
  evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);
  art::Handle< std::vector<t962::ScanInfo> > scanListHandle;
  evt.getByLabel(fScanModuleLabel,scanListHandle);
  
  art::Handle< std::vector<recob::Cluster>  > kingaListHandle;
  evt.getByLabel(fKingaModuleLabel,kingaListHandle);
  art::Handle< std::vector<recob::Cluster> > linemergerclusterListHandle;
  evt.getByLabel(fLineMergerModuleLabel,linemergerclusterListHandle);
  art::Handle< std::vector<recob::Cluster> > dbscanclusterListHandle;
  evt.getByLabel(fDbscanModuleLabel,dbscanclusterListHandle);

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

  art::PtrVector<t962::MINOS> minoslist;
  if(evt.getByLabel(fMINOSModuleLabel,minosListHandle))
  for (unsigned int i = 0; i < minosListHandle->size(); i++){
    art::Ptr<t962::MINOS> minosHolder(minosListHandle,i);
    minoslist.push_back(minosHolder);
  }

  art::PtrVector<t962::MINOSTrackMatch> trackmatchlist;
  if(evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle))
  for (unsigned int i = 0; i < trackmatchListHandle->size(); i++){
    art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,i);
    trackmatchlist.push_back(trackmatchHolder);
  }
  
   art::PtrVector<t962::ScanInfo> scanlist;
  if(evt.getByLabel(fScanModuleLabel,scanListHandle))
  for (unsigned int i = 0; i < scanListHandle->size(); i++){
    art::Ptr<t962::ScanInfo> scanHolder(scanListHandle,i);
    scanlist.push_back(scanHolder);
  }
  art::ServiceHandle<geo::Geometry> geom;  
  art::ServiceHandle<util::LArProperties> larp;

  //vertex information
  if(vertexlist.size())
  {
    double vtxxyz[3];
    vertexlist[0]->XYZ(vtxxyz);
    vtxx_reco = vtxxyz[0];
    vtxy_reco = vtxxyz[1];
    vtxz_reco = vtxxyz[2];
  }
  
  
  // 2d vertex information
  bool found2dvtx = false;
  
  for (unsigned int j = 0; j<endpointlist.size();j++){
 std::cout<<"j="<<j<<" W_VERTEX_RECO= "<<endpointlist[j]->WireNum()<<" T_VERTEX_RECO= "<<endpointlist[j]->DriftTime()<<std::endl;
          twodvtx_w_reco[j]=endpointlist[j]->WireNum();
          twodvtx_t_reco[j]=endpointlist[j]->DriftTime();
          found2dvtx=true;
    }
  
  //kingacluster information
  
    int fkingaCl_p0=0;
    int fkingaCl_p1=0;
    int fkingaCl_near_vertex_p0=0;
    int fkingaCl_near_vertex_p1=0;
     if(evt.getByLabel(fKingaModuleLabel,kingaListHandle)){
     art::PtrVector<recob::Cluster> KingaClusIn;  
     
     
     for(unsigned int ii = 0; ii < kingaListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(kingaListHandle, ii);
      KingaClusIn.push_back(cluster);
      
      
      if(cluster->View()==geo::kU){
      kingaclusters_planeNo[ii]=0;
      Start_pt_w_kingaCl[ii]=cluster->StartPos()[0];
      Start_pt_t_kingaCl[ii]=cluster->StartPos()[1];
      fkingaCl_p0++;
      //std::cout<<"p0, cluster# "<<ii<<" startPoint: "<<cluster->StartPos()[0]<<" , "<<cluster->StartPos()[1]<<std::endl;
       if(fabs(cluster->StartPos()[0]-twodvtx_w_reco[0])<6 && fabs(cluster->StartPos()[1]-twodvtx_t_reco[0])<90 ){
      fkingaCl_near_vertex_p0++;
      
        }
      
      }
       else if(cluster->View()==geo::kV){
      kingaclusters_planeNo[ii]=1;
      Start_pt_w_kingaCl[ii]=cluster->StartPos()[0];
      Start_pt_t_kingaCl[ii]=cluster->StartPos()[1];
       fkingaCl_p1++;
       //std::cout<<"p1, cluster# "<<ii<<" startPoint: "<<cluster->StartPos()[0]<<" , "<<cluster->StartPos()[1]<<std::endl;
         if( fabs(cluster->StartPos()[0]-twodvtx_w_reco[1])<6 && fabs(cluster->StartPos()[1]-twodvtx_t_reco[1])<90 ){
         fkingaCl_near_vertex_p1++;
      
        }
      }
   }
    no_kingaclusters=KingaClusIn.size();
   } //if
   
  nkingaclustersu_reco=fkingaCl_p0;
  nkingaclustersv_reco=fkingaCl_p1;
  nvertexkingaclustersu_reco=fkingaCl_near_vertex_p0;
  nvertexkingaclustersv_reco=fkingaCl_near_vertex_p1;
   

  //line merger cluster information
  
    int flinemergerCl_p0=0;
    int flinemergerCl_p1=0;
    int flinemergerCl_near_vertex_p0=0;
    int flinemergerCl_near_vertex_p1=0;
     if(evt.getByLabel(fLineMergerModuleLabel,linemergerclusterListHandle)){
     art::PtrVector<recob::Cluster> LineMergerClusIn;  
     for(unsigned int ii = 0; ii < linemergerclusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(linemergerclusterListHandle, ii);
      LineMergerClusIn.push_back(cluster);
      
      
      if(cluster->View()==geo::kU){
      linemergerclusters_planeNo[ii]=0;
      Start_pt_w_linemergerCl[ii]=cluster->StartPos()[0];
      Start_pt_t_linemergerCl[ii]=cluster->StartPos()[1];
      flinemergerCl_p0++;
       if(found2dvtx==true && fabs(cluster->StartPos()[0]-twodvtx_w_reco[0])<6 && fabs(cluster->StartPos()[1]-twodvtx_t_reco[0])<90 ){
      flinemergerCl_near_vertex_p0++;
      
        }
      
      }
       else if(cluster->View()==geo::kV){
       linemergerclusters_planeNo[ii]=1;
       Start_pt_w_linemergerCl[ii]=cluster->StartPos()[0];
       Start_pt_t_linemergerCl[ii]=cluster->StartPos()[1];
       flinemergerCl_p1++;
         if(found2dvtx==true && fabs(cluster->StartPos()[0]-twodvtx_w_reco[1])<6 && fabs(cluster->StartPos()[1]-twodvtx_t_reco[1])<90 ){
         flinemergerCl_near_vertex_p1++;
      
        }
      }
      }
       no_linemergerclusters=LineMergerClusIn.size();
      }
  
  nlinemergerclustersu_reco=flinemergerCl_p0;
  nlinemergerclustersv_reco=flinemergerCl_p1;
  nvertexlinemergerclustersu_reco=flinemergerCl_near_vertex_p0;
  nvertexlinemergerclustersv_reco=flinemergerCl_near_vertex_p1;
  
 
  //dbscan cluster information
  
    int fdbscanCl_p0=0;
    int fdbscanCl_p1=0;
    int fdbscanCl_near_vertex_p0=0;
    int fdbscanCl_near_vertex_p1=0;
     if(evt.getByLabel(fDbscanModuleLabel,dbscanclusterListHandle)){
     art::PtrVector<recob::Cluster> DbscanClusIn;  
     for(unsigned int ii = 0; ii < dbscanclusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(dbscanclusterListHandle, ii);
      DbscanClusIn.push_back(cluster);
      
      
      if(cluster->View()==geo::kU){
      fdbscanCl_p0++;
       if(found2dvtx==true && fabs(cluster->StartPos()[0]-twodvtx_w_reco[0])<6 && fabs(cluster->StartPos()[1]-twodvtx_t_reco[0])<90 ){
      fdbscanCl_near_vertex_p0++;
      
        }
      
      }
       else if(cluster->View()==geo::kV){
       fdbscanCl_p1++;
         if(found2dvtx==true && fabs(cluster->StartPos()[0]-twodvtx_w_reco[1])<6 && fabs(cluster->StartPos()[1]-twodvtx_t_reco[1])<90 ){
         fdbscanCl_near_vertex_p1++;
      
        }
      }
      }
      }
  
  ndbscanclustersu_reco=fdbscanCl_p0;
  ndbscanclustersv_reco=fdbscanCl_p1;
  nvertexdbscanclustersu_reco=fdbscanCl_near_vertex_p0;
  nvertexdbscanclustersv_reco=fdbscanCl_near_vertex_p1;
  
  
  

  //matching information  
  nmatched_reco = trackmatchlist.size();
  int ANTtrackID=-1;
  test_charge_minos=0.;
  
  for(unsigned int j = 0; j < minoslist.size(); j++)
     { 
    
     	if (!isdata)
     	{
        	if(minoslist[j]->fcharge<0.)
        	test_charge_minos=-1;
        }
       
     }
     
  for(unsigned int i = 0; i < trackmatchlist.size(); i++)
    {
   
     for(unsigned int j = 0; j < minoslist.size(); j++)
     {

     if(minoslist[j]->ftrkIndex==trackmatchlist[i]->fMINOStrackid)
     {
     ANTtrackID=trackmatchlist[i]->fArgoNeuTtrackid;
     
     if(minoslist[j]->ftrkcontained)
     trk_mom_minos = minoslist[j]->ftrkErange;
     else
     trk_mom_minos = minoslist[j]->ftrkmom;     
     
     trk_charge_minos = minoslist[j]->fcharge;
     trk_dcosx_minos = minoslist[j]->ftrkdcosx;
     trk_dcosy_minos = minoslist[j]->ftrkdcosy;
     trk_dcosz_minos = minoslist[j]->ftrkdcosz;
     trk_vtxx_minos = minoslist[j]->ftrkVtxX;
     trk_vtxy_minos = minoslist[j]->ftrkVtxY;
     trk_vtxz_minos = minoslist[j]->ftrkVtxZ;     
      }    
      
     }
    
    }

      //track information
     ntracks_reco=tracklist.size();
     double larStart[3];
     double larEnd[3];
 	 std::vector<double> trackStart;
  	 std::vector<double> trackEnd;
     trackStart.clear();
     trackEnd.clear();
     //grab information about whether a track is associated with a vertex and whether a track leaves the detector. these variables are powerful discriminators.
     int n_vertextracks=0;
     int n_endonboundarytracks=0;
      for(unsigned int i=0; i<tracklist.size();++i){
       tracklist[i]->Extent(trackStart,trackEnd); 
      
       
       all_trackstart_x_reco[i]=trackStart[0];
       all_trackstart_y_reco[i]=trackStart[1];
       all_trackstart_z_reco[i]=trackStart[2];
       
       //the below shouldn't be here!!!!!!!!!!!!!!!!!!!!!!!!!!
       // trackstart_x_reco=trackStart[0];
//        trackstart_y_reco=trackStart[1];
//        trackstart_z_reco=trackStart[2];
//        trackexit_x_reco=trackEnd[0];
//        trackexit_y_reco=trackEnd[1];
//        trackexit_z_reco=trackEnd[2];  



        if (!isdata && mclist[0]->NeutrinoSet()!=0){   
       
if(sqrt(pow(trackstart_x_reco-mclist[0]->GetNeutrino().Nu().Vx(),2)+pow(trackstart_y_reco-mclist[0]->GetNeutrino().Nu().Vy(),2)+pow(trackstart_z_reco-mclist[0]->GetNeutrino().Nu().Vz(),2))<fvertextrackWindow)
       n_vertextracks++; 
       }
      
       //kinga:
       if (isdata){        if(sqrt(pow(trackstart_x_reco-vtxx_reco,2)+pow(trackstart_y_reco-vtxy_reco,2)+pow(trackstart_z_reco-vtxz_reco,2))<fvertextrackWindow)
       n_vertextracks++; 
       }
       
       if(EndsOnBoundary(tracklist[i])) n_endonboundarytracks++;
      }
     
     nvertextracks_reco=n_vertextracks; 
     ntrackendonboundary_reco=n_endonboundarytracks;
     
     //grab information about where track started and ended and the dcos at those points
     trackStart.clear();
     trackEnd.clear();
     memset(larStart, 0, 3);
     memset(larEnd, 0, 3);
     for(unsigned int i=0; i<tracklist.size();++i){
     
      //kinga:
      if(tracklist.size()==2){ // (2 track event)
      
       tracklist[i]->Direction(larStart,larEnd);
      
       two_trackstart_dcosx_reco[i] = larStart[0];
       two_trackstart_dcosy_reco[i] = larStart[1];
       two_trackstart_dcosz_reco[i] = larStart[2];       
       two_trackexit_dcosx_reco[i] = larEnd[0];
       two_trackexit_dcosy_reco[i] = larEnd[1];
       two_trackexit_dcosz_reco[i] = larEnd[2];
      }
      
      
      if(ANTtrackID!=tracklist[i]->ID())
      continue;

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

    //scan information (data only)
     double time= -99.;
     double y_vert=-99.;
     double z_vert=-99.;
     
     for(unsigned int i=0; i<scanlist.size();++i){ 
     time=(scanlist[i]->Get_VertIndTime()+scanlist[i]->Get_VertColTime())/2.;  
     int wire_I=scanlist[i]->Get_VertIndWire();    
     int wire_C=scanlist[i]->Get_VertColWire();   
     if(wire_I>-1 && wire_C>-1)
     {
     wire_I=geom->PlaneWireToChannel(0,wire_I);
     wire_C=geom->PlaneWireToChannel(1,wire_C);
     geom->ChannelsIntersect(wire_I,wire_C,y_vert,z_vert);
         
     vtxx_scan=(time-60)*.031;
     vtxy_scan=y_vert;
     vtxz_scan=z_vert;
     neutrino_scan=scanlist[i]->Get_IsNeutrino();
     maybeneutrino_scan=scanlist[i]->Get_IsMaybeNeutrino();
     ntracks_scan=scanlist[i]->Get_Track();
     nshowers_scan=scanlist[i]->Get_NumShower(); 
     }
     }
    //mc truth information
   if (!isdata && mclist[0]->NeutrinoSet()!=0){
    nuPDG_truth = mclist[0]->GetNeutrino().Nu().PdgCode();
    ccnc_truth = mclist[0]->GetNeutrino().CCNC();
    mode_truth = mclist[0]->GetNeutrino().Mode();
    Q2_truth = mclist[0]->GetNeutrino().QSqr();
    W_truth = mclist[0]->GetNeutrino().W();
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
   
    //get true 2d vertex:
    
    
    for( unsigned int i = 0; i < mclist.size(); ++i ){

    art::Ptr<simb::MCTruth> mc(mclist[i]);

    simb::MCParticle neut(mc->GetParticle(i));

    
    fMCvertex[0] =neut.Vx();
    fMCvertex[1] =neut.Vy();
    fMCvertex[2] =neut.Vz();
   
    double presamplings=60.0;
    double drifttick=(fMCvertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198)+presamplings;
    
    twodvtx_t_truth[0]=drifttick;
    twodvtx_t_truth[1]=drifttick;
    
    
    
    //.....
     genie_no_primaries=mc->NParticles();
  
     for(int j = 0; j < mc->NParticles(); ++j){
    simb::MCParticle part(mc->GetParticle(j));
    
    std::cout<<"pdg= "<<part.PdgCode()<<" ,Process="<<part.Process()<<" StatusCode= "<<part.StatusCode()<<" mass= "<<part.Mass()<<" p= "<<part.P()<<" E= "<<part.E()<<" trackID= "<<part.TrackId()<<" ND= "<<part.NumberDaughters()<<" Mother= "<<part.Mother()<<std::endl;
    
    
    
  genie_primaries_pdg[j]=part.PdgCode();
genie_Eng[j]=part.E();
genie_Px[j]=part.Px();
 genie_Py[j]=part.Py();
genie_Pz[j]=part.Pz();
genie_P[j]=part.Px();
 genie_status_code[j]=part.StatusCode();
 genie_mass[j]=part.Mass();
 genie_trackID[j]=part.TrackId();
 genie_ND[j]=part.NumberDaughters();
genie_mother[j]=part.Mother();
    
    
    
    }
    
  }
  // now wire vertex:
    unsigned int channel2,plane2,wire2,tpc2,cs; 
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
  mf::LogWarning("ccqeanalysistreeexcp")<<e;
  
    
  //std::cout<<"fMCvertex[2]= "<<fMCvertex[2]<<" plane="<<plane<<" DetLength= "<<geom->DetLength()<<" geom->Nchannels()= "<<geom->Nchannels()<<std::endl;
  
  if(plane==0 && fMCvertex[2]<5) channel2=0;
  else if(plane==0 && fMCvertex[2]>geom->DetLength()-5) channel2=(geom->Nchannels())/2 -1;
  else if(plane==1 && fMCvertex[2]>geom->DetLength()-5) channel2=geom->Nchannels()-1;
  else if(plane==1 && fMCvertex[2]<5) channel2=(geom->Nchannels())/2 -1;

  
  }
   geom->ChannelToWire(channel2,cs,tpc2,plane2,wire2);   
   
   
   twodvtx_w_truth[plane]=wire2;
   }
   }
  
    
    
    
    
    
  }//MC
  
  
  
  
  
  
  
   //--------------------------------------------------------------//
 //        NOW I WILL GET INFO FROM GEANT4 TO FIND OUT HOW MANY 
 //        PARTICLES WE CAN REALLY SEE IN OUR DETECTOR
 //        this is needed if you want to confirm that kingaclusters 
 ///       can correctly count tracks:
 //--------------------------------------------------------------//
 
 
  if (!isdata){ 
 
 art::Handle< std::vector<sim::Particle> > geant_list;
   if(evt.getByLabel (fLArG4ModuleLabel,geant_list));
 
  art::PtrVector<sim::Particle> geant_part;
   for (unsigned int ii = 0; ii <  geant_list->size(); ++ii)
    {
      art::Ptr<sim::Particle> p(geant_list,ii);
      geant_part.push_back(p);
    } 
    std::cout<<"No of geant part= "<<geant_list->size()<<std::endl;
 std::string pri ("primary");
 int primary=0;
 //determine the number of primary particles from geant:
  
  for( unsigned int i = 0; i < geant_part.size(); ++i ){
   
    if(geant_part[i]->Process()==pri){
    primary++;
    }
   
   }
  
  no_primaries=primary;
  
 std::cout<<"Geant4 list: "<<std::endl;
 
 for( unsigned int i = 0; i < geant_part.size(); ++i ){
   std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" Process= "<<geant_part[i]->Process()<<" E= "<<geant_part[i]->E()<<" P= "<<geant_part[i]->P()<<" "<<sqrt(geant_part[i]->Px()*geant_part[i]->Px() + geant_part[i]->Py()*geant_part[i]->Py()+ geant_part[i]->Pz()*geant_part[i]->Pz())<<std::endl;
   
   if(geant_part[i]->Process()==pri){
   
   //std::cout<<"StatusCode= "<<geant_part[i]->StatusCode()<<" Mother= "<<geant_part[i]->Mother()<<std::endl;
   
    // fprimaries_pdg.push_back(geant_part[i]->PdgCode());
//     std::cout<<"geant_part[i]->E()= "<<geant_part[i]->E()<<std::endl;
//     fEng.push_back(geant_part[i]->E());
//     
   
    primaries_pdg[i]=geant_part[i]->PdgCode();
    
    Eng[i]=geant_part[i]->E();
    Px[i]=geant_part[i]->Px();
   
    Py[i]=geant_part[i]->Py();
    Pz[i]=geant_part[i]->Pz();
    
   StartPointx[i]=geant_part[i]->Vx();
   StartPointy[i]=geant_part[i]->Vy();
   StartPointz[i]=geant_part[i]->Vz();
   EndPointx[i]=geant_part[i]->EndPoint()[0];
   EndPointy[i]=geant_part[i]->EndPoint()[1];
   EndPointz[i]=geant_part[i]->EndPoint()[2];
   
   NumberDaughters[i]=geant_part[i]->NumberDaughters();
   
   
   std::cout<<"length= "<<sqrt((EndPointx[i]-StartPointx[i])*(EndPointx[i]-StartPointx[i]) + (EndPointy[i]-StartPointy[i])*(EndPointy[i]-StartPointy[i])+ (EndPointz[i]-StartPointz[i])*(EndPointz[i]-StartPointz[i]))<<std::endl;
   
 // std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" trackId= "<<geant_part[i]->TrackId()<<" mother= "<<geant_part[i]->Mother()<<" NumberDaughters()= "<<geant_part[i]->NumberDaughters()<<" process= "<<geant_part[i]->Process()<<std::endl;
     
     }
     
     }
 
 
} //if MC

//////////////////////////////////////////////////////////////////////////
/////////        HIT info:       ////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
  std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
  unsigned int p(0),w(0),t(0), channel(0),cs(0);
   int hit_no;
   
   
   std::cout<<"hits.size()= "<<hits.size()<<std::endl;
   
   no_hits=hits.size();
   
  while( itr != hits.end() ){
 // std::cout<<"working on hit# "<<itr-hits.begin()<<std::endl;
  hit_no=(int)(itr-hits.begin());
  //std::cout<<"hit_no= "<<hit_no<<std::endl;
  //std::cout<<"channel= "<<(*itr)->Channel()<<" ";
  channel=(*itr)->Channel();
  
  geom->ChannelToWire(channel,cs,t,p,w);
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
}

  //---------------------------------------------------------------- 
void t962::CCQEAnalysisTree::ResetVars(){

  run = -99999;
  event = -99999;
  isdata = -99999;
  vtxx_reco = -99999;
  vtxy_reco  = -99999;
  vtxz_reco  = -99999;
  
  ntracks_reco = -99999;
 
  nvertexclustersw_reco = -99999;
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
    
  vtxx_scan=-99999;
  vtxy_scan=-99999;
  vtxz_scan=-99999;
  neutrino_scan=-99999;
  maybeneutrino_scan=-99999;
  ntracks_scan=-99999;
  nshowers_scan=-99999;         
  nuPDG_truth = -99999;
  ccnc_truth = -99999;
  mode_truth = -99999;
  enu_truth = -99999;
  Q2_truth = -99999;
  W_truth = -99999;
  nuvtxx_truth = -99999;
  nuvtxy_truth = -99999;
  nuvtxz_truth = -99999;
  lep_mom_truth = -99999;
  lep_dcosx_truth = -99999;
  lep_dcosy_truth = -99999;
  lep_dcosz_truth = -99999;
  
  nkingaclustersu_reco=-9999;
  nkingaclustersv_reco=-9999;
  nvertexkingaclustersu_reco=-9999;
  nvertexkingaclustersv_reco=-9999;
  nlinemergerclustersu_reco=-9999;
  nlinemergerclustersv_reco=-9999;
  nvertexlinemergerclustersu_reco=-9999;
  nvertexlinemergerclustersv_reco=-9999;
  ndbscanclustersu_reco=-9999;
  ndbscanclustersv_reco=-9999;
  nvertexdbscanclustersu_reco=-9999;
  nvertexdbscanclustersv_reco=-9999;
  no_kingaclusters=-999;
  no_linemergerclusters=-999;
   ntracks_reco=-999;
   no_primaries=-999;
   genie_no_primaries=-999;
   no_hits=-999;
  
}

bool t962::CCQEAnalysisTree::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
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








