////////////////////////////////////////////////////////////////////////
// \version $Id$
// 
// \author tjyang@fnal.gov
// \author joshua.spitz@yale.edu
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TTree;

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxClusters   = 1000;  //maximum number of clusters
const int kMaxHits       = 20000; //maximum number of hits;
const int kMaxPrimaries  = 1000;  //maximum number of primary particles
const int kMaxTrackHits  = 1000;  //maximum number of hits on a track

namespace recob{
  class Track;
  class Hit;
}

///Track finding and building 
namespace t962 {
   
  class AnalysisTree : public art::EDAnalyzer {

  public:
          
    explicit AnalysisTree(fhicl::ParameterSet const& pset); 
    virtual ~AnalysisTree();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);
    
  protected:
  
    bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);
  
  private:
    
    void HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, int& trackid, double& purity, double& maxe);
    void ResetVars();

    TTree* fTree;
    //run information
    int run;
    int subrun;
    int event;
    double evttime;
    double beamtime;
    double pot;
    double taulife;
    int isdata;
    //reconstructed information    
    double vtxx_reco;
    double vtxy_reco;
    double vtxz_reco;
    int nclusu_reco;
    int nclusv_reco;
    int nclusw_reco;
    int ntracks_reco;         //number of reconstructed tracks
    int nvertextracks_reco;   //number of reconstructed tracks with start position within fvertextrackWindow cm of vertex
    int ntrackendonboundary_reco; //number of reconstructed tracks with end poistion within fboundaryWindow cm of detector boundary
    int nvertexclustersu_reco; //number of reconstructed clusters with start position within fvertexclusterWindow cm of vertex
    int nvertexclustersv_reco; //number of reconstructed clusters with start position within fvertexclusterWindow cm of vertex
    int nvertexclustersw_reco; //number of reconstructed clusters with start position within fvertexclusterWindow cm of vertex
    double trackstart_x_reco;
    double trackstart_y_reco;
    double trackstart_z_reco;
    double trackexit_x_reco;
    double trackexit_y_reco;
    double trackexit_z_reco;
    double trackstart_dcosx_reco;
    double trackstart_dcosy_reco;
    double trackstart_dcosz_reco;       
    double trackexit_dcosx_reco;
    double trackexit_dcosy_reco;
    double trackexit_dcosz_reco;   
    double enu_reco;
    int    nclupertrack_reco;  //how many times one cluster is used in tracking?
    double trkvtxx[kMaxTrack];
    double trkvtxy[kMaxTrack];
    double trkvtxz[kMaxTrack];
    double trkendx[kMaxTrack];
    double trkendy[kMaxTrack];
    double trkendz[kMaxTrack];
    double trkstartdcosx[kMaxTrack];
    double trkstartdcosy[kMaxTrack];
    double trkstartdcosz[kMaxTrack];
    double trkenddcosx[kMaxTrack];
    double trkenddcosy[kMaxTrack];
    double trkenddcosz[kMaxTrack];
    double trkke[kMaxTrack];
    double Kin_Eng_reco[kMaxTrack];
    double trkrange[kMaxTrack];
    double trkpitchc[kMaxTrack];
    int    nmaxtrkhits;
    int    ntrkhits[kMaxTrack];
    double trkdedx[kMaxTrack][kMaxTrackHits];
    double trkdqdx[kMaxTrack][kMaxTrackHits];
    double trkresrg[kMaxTrack][kMaxTrackHits];
    int    trkhitindexc[kMaxTrack][kMaxTrackHits];
    int    trkpid[kMaxTrack];
    int    trkpidndf[kMaxTrack];
    double trkpidchi2[kMaxTrack];
    double trkmissinge[kMaxTrack];
    double trkmissingeavg[kMaxTrack];
    int    trktruepdgu[kMaxTrack];
    double trktrueeffu[kMaxTrack];
    double trktruepuru[kMaxTrack];
    int    trktruepdgv[kMaxTrack];
    double trktrueeffv[kMaxTrack];
    double trktruepurv[kMaxTrack];
    //matching information
    int nmatched_reco;        //number of matched tracks
    double trk_mom_minos;
    double trk_charge_minos;
    double trk_dcosx_minos;
    double trk_dcosy_minos;
    double trk_dcosz_minos;
    double trk_vtxx_minos;
    double trk_vtxy_minos;
    double trk_vtxz_minos;
    double trk_index_minos;
    int trkcontained_minos;
    int test_charge_minos;
    double rdiff_minos;
    double thetadiff_minos;
    double muon_Kin_Eng_reco;
    int no_dead_wires_muon;
//    //scan information
//    double vtxx_scan;
//    double vtxy_scan;
//    double vtxz_scan;
//    int neutrino_scan;
//    int maybeneutrino_scan;
//    int ntracks_scan;
//    int nshowers_scan; 
    
    //mctruth information
    int parpdg;
    double parmom;
    int nuPDG_truth;
    int ccnc_truth;
    int mode_truth;
    double enu_truth;
    double Q2_truth;
    double W_truth;
    double nuvtxx_truth;
    double nuvtxy_truth;
    double nuvtxz_truth;
    double nu_dcosx_truth;
    double nu_dcosy_truth;
    double nu_dcosz_truth;
    double lep_mom_truth;
    double lep_dcosx_truth;
    double lep_dcosy_truth;
    double lep_dcosz_truth;
    int    mcevts_truth;
    double beamwgt;
    double tpx_flux;
    double tpy_flux;
    double tpz_flux;
    int    tptype_flux;
    
    float mc_index_minos;
    double mc_pdg_minos;
    double mc_px_minos;
    double mc_py_minos;
    double mc_pz_minos;
    double mc_ene_minos;
    double mc_mass_minos;
    double mc_vtxx_minos;
    double mc_vtxy_minos;
    double mc_vtxz_minos;
    int hitnuc_truth;

    //Kinga
    double twodvtx_w_reco[2];
    double twodvtx_t_reco[2];
    double twodvtx_w_truth[2];
    double twodvtx_t_truth[2];

    int nkingaclustersu_reco;
    int nkingaclustersv_reco;
    int nvertexkingaclustersu_reco;
    int nvertexkingaclustersv_reco;
    
    int nlinemergerclustersu_reco;
    int nlinemergerclustersv_reco;
    int nvertexlinemergerclustersu_reco;
    int nvertexlinemergerclustersv_reco;
    
    int ndbscanclustersu_reco;
    int ndbscanclustersv_reco;
    int nvertexdbscanclustersu_reco;
    int nvertexdbscanclustersv_reco;
    
    int no_kingaclusters;
    int no_linemergerclusters;

    int kingaclusters_planeNo[kMaxClusters];
    double Start_pt_w_kingaCl[kMaxClusters];
    double Start_pt_t_kingaCl[kMaxClusters];
    int linemergerclusters_planeNo[kMaxClusters];
    double Start_pt_w_linemergerCl[kMaxClusters];
    double Start_pt_t_linemergerCl[kMaxClusters];
    
    double two_trackstart_dcosx_reco[2];
    double two_trackstart_dcosy_reco[2];
    double two_trackstart_dcosz_reco[2];
    double two_trackexit_dcosx_reco[2];
    double two_trackexit_dcosy_reco[2];
    double two_trackexit_dcosz_reco[2];

    int no_primaries;
    int geant_list_size;
    int pdg[kMaxPrimaries];
    double Eng[kMaxPrimaries];
    double Px[kMaxPrimaries];
    double Py[kMaxPrimaries];
    double Pz[kMaxPrimaries];
    double StartPointx[kMaxPrimaries];
    double StartPointy[kMaxPrimaries];
    double StartPointz[kMaxPrimaries];
    double EndPointx[kMaxPrimaries];
    double EndPointy[kMaxPrimaries];
    double EndPointz[kMaxPrimaries];
    int NumberDaughters[kMaxPrimaries];
    int TrackId[kMaxPrimaries];
    int Mother[kMaxPrimaries];
    int process_primary[kMaxPrimaries];
    //from genie:
    int genie_no_primaries;
    int genie_primaries_pdg[kMaxPrimaries];
    double genie_Eng[kMaxPrimaries];
    double genie_Px[kMaxPrimaries];
    double genie_Py[kMaxPrimaries];
    double genie_Pz[kMaxPrimaries];
    double genie_P[kMaxPrimaries];
    int genie_status_code[kMaxPrimaries];
    double genie_mass[kMaxPrimaries];
    int genie_trackID[kMaxPrimaries];
    int genie_ND[kMaxPrimaries];
    int genie_mother[kMaxPrimaries];
       
    int no_hits;
    int hit_plane[kMaxHits];
    int hit_wire[kMaxHits];
    int hit_channel[kMaxHits];
    double hit_peakT[kMaxHits];
    double hit_charge[kMaxHits];
    double hit_ph[kMaxHits];
    double hit_etruth[kMaxHits];

    //paddle information
    double pmttime;
    int    pmt1[4];
    int    pmt2[4];
    int    pmt3[4];
    int    pmt4[4];

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel; 
    std::string fGenieGenModuleLabel;
    std::string fG4ModuleLabel;
    std::string fClusterModuleLabel;   
    std::string fKingaModuleLabel;
    std::string fLineMergerModuleLabel;
    std::string fDbscanModuleLabel;
    std::string fTrackModuleLabel;
    std::string fEndPoint2DModuleLabel;
    std::string fVertexModuleLabel;
    std::string fMINOSModuleLabel;
    std::string fTrackMatchModuleLabel;
    //std::string fScanModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    std::string fPaddlesModuleLabel;
    double fvertextrackWindow;
    double fvertexclusterWindow;
    double fboundaryWindow;

    TH1D *hBeamWeight_numu_numode;

  }; // class AnalysisTree
}

#endif 
