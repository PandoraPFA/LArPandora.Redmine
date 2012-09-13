////////////////////////////////////////////////////////////////////////
// \version $Id$
//
// \brief Create a TTree for analysis
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
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TString.h"
#include "TTimeStamp.h"

#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Principal/View.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/FindMany.h"


#include "T962/AnalysisTree/AnalysisTree.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "AnalysisBase/anabase.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/summary.h"
#include "Simulation/SimListUtils.h"
#include "MCCheater/BackTracker.h"
 
const double MPi  = 0.1396;
const double MPro = 0.938272;
const double MMu  = 0.1057;
double LifetimeCorrection(float time){

   double t = time;

   art::ServiceHandle<util::LArProperties> LArProp;
   art::ServiceHandle<util::DetectorProperties> detprop;

   double timetick = detprop->SamplingRate()*1.e-3;    //time sample in microsec
   double presamplings = detprop->TriggerOffset();

   t -= presamplings;
   time = t * timetick;  //  (in microsec)

   double tau = LArProp->ElectronLifetime();

   double correction = exp(time/tau);
   return correction;
}


//-------------------------------------------------

t962::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) : 
  
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")      ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fG4ModuleLabel            (pset.get< std::string >("G4ModuleLabel")           ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel")      ),
  fKingaModuleLabel         (pset.get< std::string >("KingaModuleLabel")        ),
  fLineMergerModuleLabel    (pset.get< std::string >("LineMergerModuleLabel")   ),
  fDbscanModuleLabel        (pset.get< std::string >("DbscanModuleLabel")       ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")        ),
  fEndPoint2DModuleLabel    (pset.get< std::string >("EndPoint2DModuleLabel")   ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")       ),
  fMINOSModuleLabel         (pset.get< std::string >("MINOSModuleLabel")        ),
  fTrackMatchModuleLabel    (pset.get< std::string >("TrackMatchModuleLabel")   ),
  //fScanModuleLabel          (pset.get< std::string >("ScanModuleLabel")         ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel")  ),
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel")   ),
  fvertextrackWindow        (pset.get< double >("vertextrackWindow")            ),
  fvertexclusterWindow      (pset.get< double >("vertexclusterWindow")          ),
  fboundaryWindow           (pset.get< double >("boundaryWindow")               )
  
{
}

//-------------------------------------------------
t962::AnalysisTree::~AnalysisTree()
{
}

void t962::AnalysisTree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("pot",&pot,"pot/D");
  fTree->Branch("isdata",&isdata,"isdata/I");
  fTree->Branch("vtxx_reco",&vtxx_reco,"vtxx_reco/D");
  fTree->Branch("vtxy_reco",&vtxy_reco,"vtxy_reco/D");
  fTree->Branch("vtxz_reco",&vtxz_reco,"vtxz_reco/D");
  fTree->Branch("nclusu_reco",&nclusu_reco,"nclusu_reco/I");
  fTree->Branch("nclusv_reco",&nclusv_reco,"nclusv_reco/I");
  fTree->Branch("nclusw_reco",&nclusw_reco,"nclusw_reco/I");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("nvertextracks_reco",&nvertextracks_reco,"nvertextracks_reco/I");
  fTree->Branch("nvertexclustersu_reco",&nvertexclustersu_reco,"nvertexclustersu_reco/I");
  fTree->Branch("nvertexclustersv_reco",&nvertexclustersv_reco,"nvertexclustersv_reco/I");
  fTree->Branch("nvertexclustersw_reco",&nvertexclustersw_reco,"nvertexclustersw_reco/I");
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
  fTree->Branch("enu_reco",&enu_reco,"enu_reco/D");
  fTree->Branch("nclupertrack_reco",&nclupertrack_reco,"nclupertrack_reco/I");
  fTree->Branch("trkvtxx",trkvtxx,"trkvtxx[ntracks_reco]/D");
  fTree->Branch("trkvtxy",trkvtxy,"trkvtxy[ntracks_reco]/D");
  fTree->Branch("trkvtxz",trkvtxz,"trkvtxz[ntracks_reco]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("trkke",trkke,"trkke[ntracks_reco]/D");
  fTree->Branch("trkrange",trkrange,"trkrange[ntracks_reco]/D");
  fTree->Branch("trkpitchc",trkpitchc,"trkpitchc[ntracks_reco]/D");
  fTree->Branch("trkpid",trkpid,"trkpid[ntracks_reco]/I");
  fTree->Branch("trkpidndf",trkpidndf,"trkpidndf[ntracks_reco]/I");
  fTree->Branch("trkpidchi2",trkpidchi2,"trkpidchi2[ntracks_reco]/D");
  fTree->Branch("trkmissinge",trkmissinge,"trkmissinge[ntracks_reco]/D");
  fTree->Branch("trkmissingeavg",trkmissingeavg,"trkmissingeavg[ntracks_reco]/D");
  fTree->Branch("trktruepdgu",trktruepdgu,"trktruepdgu[ntracks_reco]/I");
  fTree->Branch("trktrueeffu",trktrueeffu,"trktrueeffu[ntracks_reco]/D");
  fTree->Branch("trktruepuru",trktruepuru,"trktruepuru[ntracks_reco]/D");
  fTree->Branch("trktruepdgv",trktruepdgv,"trktruepdgv[ntracks_reco]/I");
  fTree->Branch("trktrueeffv",trktrueeffv,"trktrueeffv[ntracks_reco]/D");
  fTree->Branch("trktruepurv",trktruepurv,"trktruepurv[ntracks_reco]/D");  
  fTree->Branch("nmatched_reco",&nmatched_reco,"nmatched_reco/I");  
  fTree->Branch("trk_mom_minos",&trk_mom_minos,"trk_mom_minos/D");
  fTree->Branch("trk_charge_minos",&trk_charge_minos,"trk_charge_minos/D");
  fTree->Branch("trk_dcosx_minos",&trk_dcosx_minos,"trk_dcosx_minos/D");
  fTree->Branch("trk_dcosy_minos",&trk_dcosy_minos,"trk_dcosy_minos/D");
  fTree->Branch("trk_dcosz_minos",&trk_dcosz_minos,"trk_dcosz_minos/D");
  fTree->Branch("trk_vtxx_minos",&trk_vtxx_minos,"trk_vtxx_minos/D");
  fTree->Branch("trk_vtxy_minos",&trk_vtxy_minos,"trk_vtxy_minos/D");
  fTree->Branch("trk_vtxz_minos",&trk_vtxz_minos,"trk_vtxz_minos/D"); 
  fTree->Branch("mc_index_minos",&mc_index_minos,"mc_index_minos/F");
  fTree->Branch("mc_pdg_minos",&mc_pdg_minos,"mc_pdg_minos/D");
  fTree->Branch("mc_px_minos",&mc_px_minos,"mc_px_minos/D");
  fTree->Branch("mc_py_minos",&mc_py_minos,"mc_py_minos/D");
  fTree->Branch("mc_pz_minos",&mc_pz_minos,"mc_pz_minos/D");
  fTree->Branch("mc_ene_minos",&mc_ene_minos,"mc_ene_minos/D");
  fTree->Branch("mc_mass_minos",&mc_mass_minos,"mc_mass_minos/D");
  fTree->Branch("mc_vtxx_minos",&mc_vtxx_minos,"mc_vtxx_minos/D");
  fTree->Branch("mc_vtxy_minos",&mc_vtxy_minos,"mc_vtxy_minos/D");
  fTree->Branch("mc_vtxz_minos",&mc_vtxz_minos,"mc_vtxz_minos/D");   
  fTree->Branch("trkcontained_minos",&trkcontained_minos,"trkcontained_minos/I");  
  fTree->Branch("test_charge_minos",&test_charge_minos,"test_charge_minos/I"); 
  fTree->Branch("rdiff_minos",&rdiff_minos,"rdiff_minos/D");
  fTree->Branch("thetadiff_minos",&thetadiff_minos,"thetadiff_minos/D");
//  fTree->Branch("vtxx_scan", &vtxx_scan, "vtxx_scan/D");
//  fTree->Branch("vtxy_scan", &vtxy_scan, "vtxy_scan/D");
//  fTree->Branch("vtxz_scan", &vtxz_scan, "vtxz_scan/D");
//  fTree->Branch("ntracks_scan", &ntracks_scan , "ntracks_scan/I");
//  fTree->Branch("nshowers_scan", &nshowers_scan , "nshowers_scan/I");
//  fTree->Branch("neutrino_scan", &neutrino_scan , "neutrino_scan/I");
//  fTree->Branch("maybeneutrino_scan", &maybeneutrino_scan , "maybeneutrino_scan/I"); 
  fTree->Branch("parpdg", &parpdg, "parpdg/I");
  fTree->Branch("parmom", &parmom, "parmom/D");
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
  fTree->Branch("nu_dcosx_truth",&nu_dcosx_truth,"nu_dcosx_truth/D");
  fTree->Branch("nu_dcosy_truth",&nu_dcosy_truth,"nu_dcosy_truth/D");
  fTree->Branch("nu_dcosz_truth",&nu_dcosz_truth,"nu_dcosz_truth/D");
  fTree->Branch("lep_mom_truth",&lep_mom_truth,"lep_mom_truth/D");
  fTree->Branch("lep_dcosx_truth",&lep_dcosx_truth,"lep_dcosx_truth/D");
  fTree->Branch("lep_dcosy_truth",&lep_dcosy_truth,"lep_dcosy_truth/D");
  fTree->Branch("lep_dcosz_truth",&lep_dcosz_truth,"lep_dcosz_truth/D");
  fTree->Branch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
  fTree->Branch("beamwgt",&beamwgt,"beamwgt/D");

  //kinga:

  fTree->Branch("no_dead_wires_muon",&no_dead_wires_muon,"no_dead_wires_muon/I");  
  fTree->Branch("no_hits",&no_hits,"no_hits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[no_hits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[no_hits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[no_hits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[no_hits]/D");

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
  
  fTree->Branch("Kin_Eng_reco",Kin_Eng_reco,"Kin_Eng_reco[ntracks_reco]/D");
  //fTree->Branch("fTrkPitchC", &fTrkPitchC, "fTrkPitchC/D");
  fTree->Branch("muon_Kin_Eng_reco",&muon_Kin_Eng_reco,"muon_Kin_Eng_reco/D");
  
 //..................................
 // now from genie:
  fTree->Branch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
  fTree->Branch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
  fTree->Branch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/D");
  fTree->Branch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/D");
  fTree->Branch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/D");
  fTree->Branch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/D");
  fTree->Branch("genie_P",genie_P,"genie_P[genie_no_primaries]/D");
  fTree->Branch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
  fTree->Branch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/D");
  fTree->Branch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
  fTree->Branch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
  fTree->Branch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");
  

  TString filename = "numu_numode_final.root";
  const char *fROOTfile = gSystem->FindFile("${SRT_PRIVATE_CONTEXT}/T962/CCInclusiveMacro/:${SRT_PUBLIC_CONTEXT}/T962/CCInclusiveMacro/",filename);

  if (!fROOTfile) throw cet::exception("AnalysisTree") << "cannot find the root file: \n" 
						       << filename
						       << "\n bail ungracefully.";

  TFile *file = TFile::Open(fROOTfile);
  hBeamWeight_numu_numode = (TH1D*)file->Get("histdiv");

}


void t962::AnalysisTree::beginSubRun(const art::SubRun& sr)
{

  art::Handle< sumdata::POTSummary > potListHandle;
  //sr.getByLabel(fPOTModuleLabel,potListHandle);
  
  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    pot=potListHandle->totpot;
  else
    pot=0.;
  
}




void t962::AnalysisTree::analyze(const art::Event& evt)
{
  
  ResetVars();
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();
  
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();

  if (evt.isRealData()){
    isdata = 1;
  }
  else isdata = 0;
  
  art::Handle< std::vector<raw::RawDigit> > rdListHandle;
  evt.getByLabel(fDigitModuleLabel,rdListHandle);

  art::Handle< std::vector<sim::SimChannel> > scListHandle;
  evt.getByLabel(fDigitModuleLabel,scListHandle);

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Wire> > wireListHandle;
  evt.getByLabel(fCalDataModuleLabel,wireListHandle);

  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);


  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::Handle< std::vector<recob::Cluster>  > kingaListHandle;
  std::vector<art::Ptr<recob::Cluster> >kingalist;  
  if(evt.getByLabel(fKingaModuleLabel,kingaListHandle))
    art::fill_ptr_vector(kingalist,kingaListHandle);

  art::Handle< std::vector<recob::Cluster> > linemergerclusterListHandle;
  std::vector<art::Ptr<recob::Cluster> >linemergerlist;  
  if(evt.getByLabel(fLineMergerModuleLabel,linemergerclusterListHandle))
    art::fill_ptr_vector(linemergerlist, linemergerclusterListHandle);

  art::Handle< std::vector<recob::Cluster> > dbscanclusterListHandle;
  std::vector<art::Ptr<recob::Cluster> >dbscanlist;
  if(evt.getByLabel(fDbscanModuleLabel,dbscanclusterListHandle))
    art::fill_ptr_vector(dbscanlist,dbscanclusterListHandle);

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
  art::fill_ptr_vector(tracklist, trackListHandle);

  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  std::vector<art::Ptr<recob::EndPoint2D> > endpointlist;
  if (evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle))
    art::fill_ptr_vector(endpointlist, endpointListHandle);

  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  std::vector<art::Ptr<recob::Vertex> > vertexlist;
  if (evt.getByLabel(fVertexModuleLabel,vertexListHandle))
  art::fill_ptr_vector(vertexlist, vertexListHandle);

  art::Handle< std::vector<t962::MINOS> > minosListHandle;
  std::vector<art::Ptr<t962::MINOS> > minoslist;
  if (evt.getByLabel(fMINOSModuleLabel,minosListHandle))
    art::fill_ptr_vector(minoslist, minosListHandle);

//  art::Handle< std::vector<t962::ScanInfo> > scanListHandle;
//  std::vector<art::Ptr<t962::ScanInfo> > scanlist;
//  if (evt.getByLabel(fScanModuleLabel,scanListHandle))
//    art::fill_ptr_vector(scanlist, scanListHandle);


  art::ServiceHandle<geo::Geometry> geom;  
  art::ServiceHandle<cheat::BackTracker> bt;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> LArProp;

  // Electronic calibration factor to convert from ADC to electrons
  double fElectronsToADC = detprop->ElectronsToADC();



  //vertex information
  if(vertexlist.size())
  {
    double vtxxyz[3];
    vertexlist[0]->XYZ(vtxxyz);
    vtxx_reco = vtxxyz[0];
    vtxy_reco = vtxxyz[1];
    vtxz_reco = vtxxyz[2];
  }
  //cluster information
  nclusu_reco = 0;
  nclusv_reco = 0;
  nclusw_reco = 0;

  int nplanes = geom->Nplanes();
  std::vector<int> Cls[nplanes];

  for(size_t i=0; i<clusterlist.size();++i){
    
    switch(clusterlist[i]->View()){
    case geo::kU :
      nclusu_reco ++;
      Cls[0].push_back(i);
      break;
    case geo::kV :
      nclusv_reco ++;
      Cls[1].push_back(i);
      break;
    case geo::kW :
      nclusw_reco ++;
      Cls[2].push_back(i);
      break;
    default :
      break;
    }
  }


  for (int i = 0; i<nplanes; i++){
    int n_vertexclusters = 0;
    if (Cls[i].size()>0){
      int vtx2d_w = -99999;
      double vtx2d_t = -99999;
      bool found2dvtx = false;
      //find 2d vertex
      for (unsigned int j = 0; j<endpointlist.size();j++){
	if (endpointlist[j]->View() == clusterlist[Cls[i][0]]->View()){
	  vtx2d_w = endpointlist[j]->WireNum();
	  vtx2d_t = endpointlist[j]->DriftTime();
	  found2dvtx = true;
	  break;
	}
      }
      if (found2dvtx){
	for (unsigned j = 0; j<Cls[i].size(); j++){
	  double w = clusterlist[Cls[i][j]]->StartPos()[0];
	  double t = clusterlist[Cls[i][j]]->StartPos()[1];
	  double dtdw = clusterlist[Cls[i][j]]->dTdW();
	  double t_vtx = t+dtdw*(vtx2d_w-w);
	  double dis = TMath::Abs(vtx2d_t-t_vtx);
	  if (dis<fvertexclusterWindow){
	    n_vertexclusters++;
	  }
	}
      }
    }
    if (i==0){
      nvertexclustersu_reco = n_vertexclusters;
    }
    else if (i==1){
      nvertexclustersv_reco = n_vertexclusters;
    }
    else if (i==2){
      nvertexclustersw_reco = n_vertexclusters;
    }
  }

  
  test_charge_minos=0.;
  for(unsigned int j = 0; j < minoslist.size(); j++)
    { 
      if (!isdata)
     	{
	  if(minoslist[j]->fcharge<0.)
	    test_charge_minos=-1;
        }
    }
  
  //track information
  ntracks_reco=tracklist.size();
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  //grab information about whether a track is associated with a vertex and whether a track leaves the detector. these variables are powerful discriminators.
  int n_vertextracks=0;
  int n_endonboundarytracks=0;
  art::FindOne<anab::Calorimetry> focal(trackListHandle, evt, fCalorimetryModuleLabel);
  art::FindOne<anab::ParticleID>  fopid(trackListHandle, evt, fParticleIDModuleLabel);
  art::FindManyP<recob::Hit>      fmht(trackListHandle, evt, fTrackModuleLabel);
  art::FindMany<recob::Track>     fmtk(hitListHandle, evt, fTrackModuleLabel);
  art::FindMany<recob::Track>     fmtkcl(clusterListHandle, evt, fTrackModuleLabel);
  
  nclupertrack_reco = 0;
  for(size_t i=0; i<clusterlist.size();++i){
    if (int(fmtkcl.at(i).size())>nclupertrack_reco ){
      nclupertrack_reco = fmtkcl.at(i).size();
    }
  }
  for(unsigned int i=0; i<tracklist.size();++i){
    trackStart.clear();
    trackEnd.clear();
    memset(larStart, 0, 3);
    memset(larEnd, 0, 3);
    tracklist[i]->Extent(trackStart,trackEnd); 
    tracklist[i]->Direction(larStart,larEnd);

    //kinga:
    if(tracklist.size()==2){ // (2 track event)
      two_trackstart_dcosx_reco[i] = larStart[0];
      two_trackstart_dcosy_reco[i] = larStart[1];
      two_trackstart_dcosz_reco[i] = larStart[2];       
      two_trackexit_dcosx_reco[i] = larEnd[0];
      two_trackexit_dcosy_reco[i] = larEnd[1];
      two_trackexit_dcosz_reco[i] = larEnd[2];
    }
           
    trackstart_x_reco=trackStart[0];
    trackstart_y_reco=trackStart[1];
    trackstart_z_reco=trackStart[2];
    trackexit_x_reco=trackEnd[0];
    trackexit_y_reco=trackEnd[1];
    trackexit_z_reco=trackEnd[2];  
    if (!isdata&&mclist[0]->NeutrinoSet()){
      if(sqrt(pow(trackstart_x_reco-mclist[0]->GetNeutrino().Nu().Vx(),2)+pow(trackstart_y_reco-mclist[0]->GetNeutrino().Nu().Vy(),2)+pow(trackstart_z_reco-mclist[0]->GetNeutrino().Nu().Vz(),2))<fvertextrackWindow)
	n_vertextracks++; 
    }
    if(EndsOnBoundary(tracklist[i])) n_endonboundarytracks++;
    trkvtxx[i]        = trackStart[0];
    trkvtxy[i]        = trackStart[1];
    trkvtxz[i]        = trackStart[2];
    trkendx[i]        = trackEnd[0];
    trkendy[i]        = trackEnd[1];
    trkendz[i]        = trackEnd[2];
    trkstartdcosx[i]  = larStart[0];
    trkstartdcosy[i]  = larStart[1];
    trkstartdcosz[i]  = larStart[2];
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    if (focal.at(i).isValid()){
      trkke[i]           = focal.at(i).ref().KineticEnergy();
      trkrange[i]        = focal.at(i).ref().Range();
      no_dead_wires_muon = focal.at(i).ref().DeadWireResRC().size();
      Kin_Eng_reco[i]    = focal.at(i).ref().KineticEnergy()+no_dead_wires_muon*focal.at(i).ref().TrkPitchC()*2.1;
    }
    if (fopid.at(i).isValid()){
      trkpid[i]         = fopid.at(i).ref().Pdg();
      trkpidndf[i]      = fopid.at(i).ref().Ndf();
      trkpidchi2[i]     = fopid.at(i).ref().MinChi2();
      trkmissinge[i]    = fopid.at(i).ref().MissingE();
      trkmissingeavg[i] = fopid.at(i).ref().MissingEavg();
    }
    if (!isdata){
      // get the hits in each view
      std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(i);
      std::vector< art::Ptr<recob::Hit> > hitsU;
      std::vector< art::Ptr<recob::Hit> > hitsV;
      for(size_t ah = 0; ah < allHits.size(); ++ah){
	if     (allHits[ah]->View() == geo::kU) hitsU.push_back(allHits[ah]);
	else if(allHits[ah]->View() == geo::kV) hitsV.push_back(allHits[ah]);
      }
      int trkid;
      double purity;
      double maxe;
      HitsPurity(hitsU, trkid, purity, maxe);
      trktruepuru[i] = purity;
      if (trkid>-1){
	const sim::Particle *particle = bt->TrackIDToParticle(trkid);
	const std::vector<sim::IDE> vide = bt->TrackIDToSimIDE(trkid);
	double tote = 0;
	for (size_t iide = 0; iide<vide.size(); ++iide){
	  tote += vide[iide].energy;
	}
	trktruepdgu[i] = particle->PdgCode();
	trktrueeffu[i] = maxe/(tote/2); //I believe tote include both induction and collection energies
      }
      HitsPurity(hitsV, trkid, purity, maxe);
      trktruepurv[i] = purity;
      if (trkid>-1){
	const sim::Particle *particle = bt->TrackIDToParticle(trkid);
	const std::vector<sim::IDE> vide = bt->TrackIDToSimIDE(trkid);
	double tote = 0;
	for (size_t iide = 0; iide<vide.size(); ++iide){
	  tote += vide[iide].energy;
	}
	trktruepdgv[i] = particle->PdgCode();
	trktrueeffv[i] = maxe/(tote/2);
      }
    }
  }
  nvertextracks_reco=n_vertextracks; 
  ntrackendonboundary_reco=n_endonboundarytracks;
  

  //From CCQEAnalysisTree
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
  for(unsigned int ii = 0; ii < kingalist.size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster = kingalist[ii];
	
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
  no_kingaclusters=kingalist.size();
  
  nkingaclustersu_reco=fkingaCl_p0;
  nkingaclustersv_reco=fkingaCl_p1;
  nvertexkingaclustersu_reco=fkingaCl_near_vertex_p0;
  nvertexkingaclustersv_reco=fkingaCl_near_vertex_p1;
   
  
  //line merger cluster information
  
  int flinemergerCl_p0=0;
  int flinemergerCl_p1=0;
  int flinemergerCl_near_vertex_p0=0;
  int flinemergerCl_near_vertex_p1=0;

  for(unsigned int ii = 0; ii < linemergerlist.size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster = linemergerlist[ii];
      
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
  no_linemergerclusters=linemergerlist.size();
  
  nlinemergerclustersu_reco=flinemergerCl_p0;
  nlinemergerclustersv_reco=flinemergerCl_p1;
  nvertexlinemergerclustersu_reco=flinemergerCl_near_vertex_p0;
  nvertexlinemergerclustersv_reco=flinemergerCl_near_vertex_p1;
  
  
  //dbscan cluster information
  
  int fdbscanCl_p0=0;
  int fdbscanCl_p1=0;
  int fdbscanCl_near_vertex_p0=0;
  int fdbscanCl_near_vertex_p1=0;

  for(unsigned int ii = 0; ii < dbscanlist.size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster = dbscanlist[ii];
      
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
  
  ndbscanclustersu_reco=fdbscanCl_p0;
  ndbscanclustersv_reco=fdbscanCl_p1;
  nvertexdbscanclustersu_reco=fdbscanCl_near_vertex_p0;
  nvertexdbscanclustersv_reco=fdbscanCl_near_vertex_p1;

//////////////////////////////////////////////////////////////////////////
/////////        HIT info:       ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

  std::vector< art::Ptr<recob::Hit> >::iterator itr = hitlist.begin();
  unsigned int p(0),w(0),t(0), cs(0);
  int hit_no;
  unsigned int channel=0;

  std::cout<<"hitlist.size()= "<<hitlist.size()<<std::endl;
   
  no_hits=hitlist.size();
   
  while( itr != hitlist.end() ){
    // std::cout<<"working on hit# "<<itr-hitlist.begin()<<std::endl;
    hit_no=(int)(itr-hitlist.begin());
    //std::cout<<"hit_no= "<<hit_no<<std::endl;
    //std::cout<<"channel= "<<(*itr)->Channel()<<" ";
    channel=(*itr)->Channel();
    
    geom->ChannelToWire(channel,cs,t,p,w);
    hit_channel[hit_no]= channel;
    hit_plane[hit_no]=p;
    hit_wire[hit_no]=w;
    hit_peakT[hit_no]=(*itr)->PeakTime();
    hit_charge[hit_no]=(*itr)->Charge();
    
    std::cout<<"p= "<<p<<" ,w= "<<w<<" ,time= "<<(*itr)->PeakTime()<<" Charge= "<<(*itr)->Charge()<<" charge(true)= "<<(*itr)->Charge(true)<<std::endl;
    
    //  if((*itr)->View()==geo::kU){
    //     std::cout<<"plane 0"<<std::endl;}
    //     else if((*itr)->View()==geo::kV){
    //     std::cout<<"plane 1"<<std::endl;}
    itr++;
  }
  
//  //scan information (data only)
//  double time= -99.;
//  double y_vert=-99.;
//  double z_vert=-99.;
//  
//  for(unsigned int i=0; i<scanlist.size();++i){ 
//    time=(scanlist[i]->Get_VertIndTime()+scanlist[i]->Get_VertColTime())/2.;  
//    int wire_I=scanlist[i]->Get_VertIndWire();    
//    int wire_C=scanlist[i]->Get_VertColWire();   
//    if(wire_I>-1 && wire_C>-1)
//      {
//	wire_I=geom->PlaneWireToChannel(0,wire_I);
//	wire_C=geom->PlaneWireToChannel(1,wire_C);
//	geom->ChannelsIntersect(wire_I,wire_C,y_vert,z_vert);
//	
//	vtxx_scan=(time-60)*.031;
//	vtxy_scan=y_vert;
//	vtxz_scan=z_vert;
//	neutrino_scan=scanlist[i]->Get_IsNeutrino();
//	maybeneutrino_scan=scanlist[i]->Get_IsMaybeNeutrino();
//	ntracks_scan=scanlist[i]->Get_Track();
//	nshowers_scan=scanlist[i]->Get_NumShower(); 
//      }
//  }

  //mc truth information
  if (!isdata){
    //save single particle information
    sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);
    for ( sim::ParticleList::const_iterator ipar = plist.begin(); ipar!=plist.end(); ++ipar){
      sim::Particle *particle = ipar->second;
      parpdg = particle->PdgCode();
      parmom = particle->Momentum().P();
      break;
    }
    //save neutrino interaction information
    mcevts_truth = mclist.size();
    if (mclist[0]->NeutrinoSet()){
      //find the true neutrino corresponding to the reconstructed event
      std::map<art::Ptr<simb::MCTruth>,double> mctruthemap;
      for (size_t i = 0; i<hitlist.size(); i++){
	if (hitlist[i]->View() == geo::kV){//collection view
	  std::vector<cheat::TrackIDE> eveIDs = bt->HitToEveID(hitlist[i]);
	  for (size_t e = 0; e<eveIDs.size(); e++){
	    art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(eveIDs[e].trackID);
	    mctruthemap[mctruth]+=eveIDs[e].energy;
	  }
	}
      }
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      double maxenergy = -1;
      for (std::map<art::Ptr<simb::MCTruth>,double>::iterator ii=mctruthemap.begin(); ii!=mctruthemap.end(); ++ii){
	if ((ii->second)>maxenergy){
	  maxenergy = ii->second;
	  mctruth = ii->first;
	}
      }

      nuPDG_truth = mctruth->GetNeutrino().Nu().PdgCode();
      ccnc_truth = mctruth->GetNeutrino().CCNC();
      mode_truth = mctruth->GetNeutrino().Mode();
      Q2_truth = mctruth->GetNeutrino().QSqr();
      W_truth = mctruth->GetNeutrino().W();
      hitnuc_truth = mctruth->GetNeutrino().HitNuc();
      enu_truth = mctruth->GetNeutrino().Nu().E();
      nuvtxx_truth = mctruth->GetNeutrino().Nu().Vx();
      nuvtxy_truth = mctruth->GetNeutrino().Nu().Vy();
      nuvtxz_truth = mctruth->GetNeutrino().Nu().Vz();
      if (mctruth->GetNeutrino().Nu().P()){
	nu_dcosx_truth = mctruth->GetNeutrino().Nu().Px()/mctruth->GetNeutrino().Nu().P();
	nu_dcosy_truth = mctruth->GetNeutrino().Nu().Py()/mctruth->GetNeutrino().Nu().P();
	nu_dcosz_truth = mctruth->GetNeutrino().Nu().Pz()/mctruth->GetNeutrino().Nu().P();
      }
      lep_mom_truth = mctruth->GetNeutrino().Lepton().P();
      if (mctruth->GetNeutrino().Lepton().P()){
	lep_dcosx_truth = mctruth->GetNeutrino().Lepton().Px()/mctruth->GetNeutrino().Lepton().P();
	lep_dcosy_truth = mctruth->GetNeutrino().Lepton().Py()/mctruth->GetNeutrino().Lepton().P();
	lep_dcosz_truth = mctruth->GetNeutrino().Lepton().Pz()/mctruth->GetNeutrino().Lepton().P();
      }

      //Kinga
      double presamplings=60.0;
      double drifttick=(nuvtxx_truth/LArProp->DriftVelocity(LArProp->Efield(),LArProp->Temperature()))*(1./.198)+presamplings;
      
      twodvtx_t_truth[0]=drifttick;
      twodvtx_t_truth[1]=drifttick;

      genie_no_primaries=mctruth->NParticles();
  
      for(int j = 0; j < mctruth->NParticles(); ++j){
	simb::MCParticle part(mctruth->GetParticle(j));
    
	//std::cout<<"pdg= "<<part.PdgCode()<<" ,Process="<<part.Process()<<" StatusCode= "<<part.StatusCode()<<" mass= "<<part.Mass()<<" p= "<<part.P()<<" E= "<<part.E()<<" trackID= "<<part.TrackId()<<" ND= "<<part.NumberDaughters()<<" Mother= "<<part.Mother()<<std::endl;
    
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
      
      double fMCvertex[3];
      fMCvertex[0] = nuvtxx_truth;
      fMCvertex[1] = nuvtxy_truth;
      fMCvertex[2] = nuvtxz_truth;
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

      //--------------------------------------------------------------//
      //        NOW I WILL GET INFO FROM GEANT4 TO FIND OUT HOW MANY 
      //        PARTICLES WE CAN REALLY SEE IN OUR DETECTOR
      //        this is needed if you want to confirm that kingaclusters 
      ///       can correctly count tracks:
      //--------------------------------------------------------------//
  
      art::Handle< std::vector<sim::Particle> > geant_list;
      std::vector<art::Ptr<sim::Particle> > geant_part;
      if(evt.getByLabel (fLArG4ModuleLabel,geant_list))
	art::fill_ptr_vector(geant_part,geant_list);
 
      std::cout<<"No of geant part= "<<geant_list->size()<<std::endl;
      std::string pri("primary");
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
 

      //beam weight, only available for nu-data
      beamwgt = 1.;
      if (nuPDG_truth == 14){
	int bin = hBeamWeight_numu_numode->FindBin(enu_truth);
	if (bin>=1&&bin<=hBeamWeight_numu_numode->GetNbinsX()) 
	  beamwgt = hBeamWeight_numu_numode->GetBinContent(bin);
      }
    }
  }
  //minos matching information
  try{
  //find matched MINOS information for each track
    nmatched_reco = 0;
    art::FindOne<t962::MINOS> fomatch(trackListHandle, evt, fTrackMatchModuleLabel);
    
    //grab information about where track started and ended and the dcos at those points
    trackStart.clear();
    trackEnd.clear();
    memset(larStart, 0, 3);
    memset(larEnd, 0, 3);
    int bestmatch = -1;
    double totaldiffmin = 1e10;
    for(size_t i=0; i<tracklist.size();++i){
      
      if(!fomatch.at(i).isValid()) continue;//No matching MINOS track
      ++nmatched_reco;

      tracklist[i]->Direction(larStart,larEnd);
      tracklist[i]->Extent(trackStart,trackEnd);  

      double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane 
                                           //(this minus number is the one we measured with Mitch)
      
      double x_offset=117.4; // previously 116.9;
      double y_offset=19.3; // previously  20.28;

      double dz = D - trackEnd[2]+(100.0 * fomatch.at(i).ref().ftrkVtxZ);//z-difference between end of T962 track and
                                                                         //begin of MINOS track...in centimeters

      double l = dz/(larEnd[2]);//3-d distance between end of T962 track and begin of MINOS track

      double x_pred = l*larEnd[0]+trackEnd[0];//predicted x-pos. of T962 track at z-position equal to
                                              //start of MINOS track
      double y_pred = l*larEnd[1]+trackEnd[1];//predicted y-pos. of T962 track at z-position equal to
                                              //start of MINOS track

      double dx = 100.0*fomatch.at(i).ref().ftrkVtxX - x_offset - x_pred;
      double dy = 100.0*fomatch.at(i).ref().ftrkVtxY + y_offset - y_pred;

      double rdiff = sqrt(dx*dx + dy*dy);
      double thetadiff = TMath::ACos((larEnd[0]*fomatch.at(i).ref().ftrkdcosx)+(larEnd[1]*fomatch.at(i).ref().ftrkdcosy)+(larEnd[2]*fomatch.at(i).ref().ftrkdcosz));
      //totaldiff is a measure of the agreement between the ArgoNeuT projected track and the MINOS track based on radial distance and angle. totaldiff= rdiff/cos(theta)  
      double totaldiff=fabs(rdiff/((larEnd[0]*fomatch.at(i).ref().ftrkdcosx)+(larEnd[1]*fomatch.at(i).ref().ftrkdcosy)+(larEnd[2]*fomatch.at(i).ref().ftrkdcosz)));
      if (totaldiff<totaldiffmin){
	totaldiffmin = totaldiff;
	bestmatch = i;
	rdiff_minos = rdiff;
	thetadiff_minos = thetadiff;
      }
    }

    if (bestmatch>-1){
      trackStart.clear();
      trackEnd.clear();
      memset(larStart, 0, 3);
      memset(larEnd, 0, 3);
      
      tracklist[bestmatch]->Direction(larStart,larEnd);
      tracklist[bestmatch]->Extent(trackStart,trackEnd);  
      
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
      
      if(fomatch.at(bestmatch).ref().ftrkcontained)
	trk_mom_minos = fomatch.at(bestmatch).ref().ftrkErange;
      else
	trk_mom_minos = fomatch.at(bestmatch).ref().ftrkmom;     
      
      trkcontained_minos = fomatch.at(bestmatch).ref().ftrkcontained; 
      trk_charge_minos = fomatch.at(bestmatch).ref().fcharge;
      trk_dcosx_minos = fomatch.at(bestmatch).ref().ftrkdcosx;
      trk_dcosy_minos = fomatch.at(bestmatch).ref().ftrkdcosy;
      trk_dcosz_minos = fomatch.at(bestmatch).ref().ftrkdcosz;
      trk_vtxx_minos = fomatch.at(bestmatch).ref().ftrkVtxX;
      trk_vtxy_minos = fomatch.at(bestmatch).ref().ftrkVtxY;
      trk_vtxz_minos = fomatch.at(bestmatch).ref().ftrkVtxZ;        

      muon_Kin_Eng_reco = Kin_Eng_reco[bestmatch];

      if (!isdata){       
	mc_index_minos = fomatch.at(bestmatch).ref().fmcIndex;
	mc_pdg_minos = fomatch.at(bestmatch).ref().fmcPDG;
	mc_px_minos = fomatch.at(bestmatch).ref().fmcPx;
	mc_py_minos = fomatch.at(bestmatch).ref().fmcPy;
	mc_pz_minos = fomatch.at(bestmatch).ref().fmcPz;
	mc_ene_minos = fomatch.at(bestmatch).ref().fmcEne;
	mc_mass_minos = fomatch.at(bestmatch).ref().fmcMass;
	mc_vtxx_minos = fomatch.at(bestmatch).ref().fmcVtxX;
	mc_vtxy_minos = fomatch.at(bestmatch).ref().fmcVtxY;
	mc_vtxz_minos = fomatch.at(bestmatch).ref().fmcVtxZ;
      }

      //Calculate event energy
      enu_reco = 0;
      //sum the energy of hits not associated with tracks
      for(size_t i=0; i<hitlist.size();++i){
	if (hitlist[i]->View() == geo::kV){//collection view
	  //std::cout<<i<<" "<<fmtk.at(i).size()<<std::endl;
	  if (fmtk.at(i).size()==0){//hit not on any tracks
	    double time   = hitlist[i]->PeakTime();
	    double MIPs   = hitlist[i]->Charge(true);   // in ADC
	    double WirePitch = geom->WirePitch(0,1,1);
	    double dQdx = MIPs/WirePitch;
	    double dQdx_e = dQdx/fElectronsToADC;  // Conversion from ADC/cm to e/cm
	    
	    dQdx_e *= LifetimeCorrection(time);   // Lifetime Correction (dQdx_e in e/cm)
	    
	    double dEdx = LArProp->BirksCorrection(dQdx_e);   // Correction for charge quenching (Recombination) dEdx in MeV/cm
	    enu_reco += dEdx*WirePitch;
	  }
	}
      }
      enu_reco /= 1000;
      enu_reco = enu_reco + MPi; //assuming pion mass
      
      //sum track energy except the muon
      for(size_t i=0; i<tracklist.size();++i){
	if (int(i)!=bestmatch){ //not the matched track
	  double mass = MPi;  //default is pion mass
	  if (trkpid[i] == 2212
	      &&trkendx[i]>3&&trkendx[i]<44
	      &&trkendy[i]>-16&&trkendy[i]<16
	      &&trkendz[i]>6&&trkendz[i]<86) mass = MPro; //proton 
	  if (trkke[i]>0){
	    enu_reco += trkke[i]/1000+mass;
	  }
	}
      }
      
      //add muon energy
      enu_reco += sqrt(pow(trk_mom_minos,2)+pow(MMu,2));
    }//if bestmatch>-1
  }   
  catch( cet::exception &e){
    mf::LogWarning("AnalysisTree") << "caught exception " << e;
  }

  fTree->Fill();
}

void t962::AnalysisTree::HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, int& trackid, double& purity, double& maxe){

  trackid = -1;
  purity = -1;
  
  art::ServiceHandle<cheat::BackTracker> bt;
  
  std::map<int,double> trkide;

  for(size_t h = 0; h < hits.size(); ++h){

    art::Ptr<recob::Hit> hit = hits[h];
    std::vector<sim::IDE> ides;
    //bt->HitToSimIDEs(hit,ides);
    std::vector<cheat::TrackIDE> eveIDs = bt->HitToEveID(hit);

    for(size_t e = 0; e < eveIDs.size(); ++e){
      //std::cout<<h<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
      trkide[eveIDs[e].trackID] += eveIDs[e].energy;
    }
  }

  maxe = -1;
  double tote = 0;
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      trackid = ii->first;
    }
  }
  
  
  if (tote>0){
    purity = maxe/tote;
  }
}

  //---------------------------------------------------------------- 
void t962::AnalysisTree::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  isdata = -99999;
  vtxx_reco = -99999;
  vtxy_reco  = -99999;
  vtxz_reco  = -99999;
  nclusu_reco  = -99999;
  nclusv_reco  = -99999;
  nclusw_reco  = -99999;
  ntracks_reco = -99999;
  nvertextracks_reco = -99999;
  nvertexclustersu_reco = -99999;
  nvertexclustersv_reco = -99999;
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
  enu_reco = -99999;
  nclupertrack_reco = -99999;
  nmatched_reco = -99999;
  trk_mom_minos = -99999;
  trk_charge_minos = -99999;
  trk_dcosx_minos = -99999;
  trk_dcosy_minos = -99999;
  trk_dcosz_minos = -99999;  
  trk_vtxx_minos = -99999;
  trk_vtxy_minos = -99999;
  trk_vtxz_minos = -99999;
  muon_Kin_Eng_reco = -99999;
  no_dead_wires_muon = -99999;
  mc_index_minos=-99999;
  mc_pdg_minos=-99999;
  mc_px_minos=-99999;
  mc_py_minos=-99999;
  mc_pz_minos=-99999;
  mc_ene_minos=-99999;
  mc_mass_minos=-99999;
  mc_vtxx_minos=-99999;
  mc_vtxy_minos=-99999;
  mc_vtxz_minos=-99999;
  test_charge_minos=-99999;
  trkcontained_minos=-99999;
  rdiff_minos = 99999;
  thetadiff_minos = 99999;
//  vtxx_scan=-99999;
//  vtxy_scan=-99999;
//  vtxz_scan=-99999;
//  neutrino_scan=-99999;
//  maybeneutrino_scan=-99999;
//  ntracks_scan=-99999;
//  nshowers_scan=-99999;         
  parpdg = -99999;
  parmom = -99999;
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
  nu_dcosx_truth = -99999;
  nu_dcosy_truth = -99999;
  nu_dcosz_truth = -99999;
  lep_mom_truth = -99999;
  lep_dcosx_truth = -99999;
  lep_dcosy_truth = -99999;
  lep_dcosz_truth = -99999;
  mcevts_truth = -99999;
  beamwgt = -99999;
  for (int i = 0; i < kMaxTrack; i++){
    trkvtxx[i] = -99999;
    trkvtxy[i] = -99999;
    trkvtxz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    trkke[i] = -99999;
    Kin_Eng_reco[i] = -99999;
    trkrange[i] = -99999;
    trkpitchc[i] = -99999;
    trkpid[i] = -99999;
    trkpidndf[i] = -99999;
    trkpidchi2[i] = -99999;
    trkmissinge[i] = -99999;
    trkmissingeavg[i] = -99999;
    trktruepdgu[i] = -99999;
    trktrueeffu[i] = -99999;
    trktruepuru[i] = -99999;
    trktruepdgv[i] = -99999;
    trktrueeffv[i] = -99999;
    trktruepurv[i] = -99999;
  }    
  
  no_primaries = -99999;
  for (int i = 0; i<kMaxPrimaries; ++i){
    primaries_pdg[i] = -99999;
    Eng[i] = -99999;
    Px[i] = -99999;
    Py[i] = -99999;
    Pz[i] = -99999;
    StartPointx[i] = -99999;
    StartPointy[i] = -99999;
    StartPointz[i] = -99999;
    EndPointx[i] = -99999;
    EndPointy[i] = -99999;
    EndPointz[i] = -99999;
    NumberDaughters[i] = -99999;
    genie_primaries_pdg[i] = -99999;
    genie_Eng[i] = -99999;
    genie_Px[i] = -99999;
    genie_Py[i] = -99999;
    genie_Pz[i] = -99999;
    genie_P[i] = -99999;
    genie_status_code[i] = -99999;
    genie_mass[i] = -99999;
    genie_trackID[i] = -99999;
    genie_ND[i] = -99999;
    genie_mother[i] = -99999;
  }    
  
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -99999;
    hit_wire[i] = -99999;
    hit_channel[i] = -99999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
  }

  for (int i = 0; i<2; ++i){
    twodvtx_w_reco[i] = -99999;
    twodvtx_t_reco[i] = -99999;
    twodvtx_w_truth[i] = -99999;
    twodvtx_t_truth[i] = -99999;
  } 

  nkingaclustersu_reco = -99999;
  nkingaclustersv_reco = -99999;
  nvertexkingaclustersu_reco = -99999;
  nvertexkingaclustersv_reco = -99999;
  
  nlinemergerclustersu_reco = -99999;
  nlinemergerclustersv_reco = -99999;
  nvertexlinemergerclustersu_reco = -99999;
  nvertexlinemergerclustersv_reco = -99999;
  
  ndbscanclustersu_reco = -99999;
  ndbscanclustersv_reco = -99999;
  nvertexdbscanclustersu_reco = -99999;
  nvertexdbscanclustersv_reco = -99999;
  
  no_kingaclusters = -99999;
  no_linemergerclusters = -99999;   

  for (int i = 0; i<kMaxClusters; ++i){
   kingaclusters_planeNo[i] = -99999;
   Start_pt_w_kingaCl[i] = -99999;
   Start_pt_t_kingaCl[i] = -99999;
   linemergerclusters_planeNo[i] = -99999;
   Start_pt_w_linemergerCl[i] = -99999;
   Start_pt_t_linemergerCl[i] = -99999;
  }
  
  for (int i = 0; i<2; ++i){
    two_trackstart_dcosx_reco[i] = -99999;
    two_trackstart_dcosy_reco[i] = -99999;
    two_trackstart_dcosz_reco[i] = -99999;
    two_trackexit_dcosx_reco[i] = -99999;
    two_trackexit_dcosy_reco[i] = -99999;
    two_trackexit_dcosz_reco[i] = -99999;
  }
 
}

bool t962::AnalysisTree::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
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








