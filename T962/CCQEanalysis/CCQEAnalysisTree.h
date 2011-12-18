////////////////////////////////////////////////////////////////////////
//
// 
// \author tjyang@fnal.gov
// kinga.partyka@yale.edu
// 
////////////////////////////////////////////////////////////////////////
#ifndef CCQEAnalysisTree_H
#define CCQEAnalysisTree_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "RecoBase/Track.h"

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TTree;

///Track finding and building 
namespace t962 {

   
  class CCQEAnalysisTree : public art::EDAnalyzer {

  public:
          
    explicit CCQEAnalysisTree(fhicl::ParameterSet const& pset); 
    virtual ~CCQEAnalysisTree();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);

    
    
    
  protected:
  
  bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);
  
  private:
    
    void ResetVars();
    double fMCvertex [3];
    
    TTree* fTree;
    //run information
    int run;
    int event;
    double pot;
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
    
    double *twodvtx_w_reco;
    double *twodvtx_t_reco;
    double *twodvtx_w_truth;
    double *twodvtx_t_truth;
    double *fStart_pt_w_kingaCl;
    double *fStart_pt_t_kingaCl;
    
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
    
    int *kingaclusters_planeNo;
    double *Start_pt_w_kingaCl;
    double *Start_pt_t_kingaCl;
    int *linemergerclusters_planeNo;
    double *Start_pt_w_linemergerCl;
    double *Start_pt_t_linemergerCl;
    
    double *two_trackstart_dcosx_reco;
    double *two_trackstart_dcosy_reco;
	double *two_trackstart_dcosz_reco;
	double *two_trackexit_dcosx_reco;
	double *two_trackexit_dcosy_reco;
	double *two_trackexit_dcosz_reco;
	double *all_trackstart_x_reco;
    double *all_trackstart_y_reco;
    double *all_trackstart_z_reco;
              
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
    int test_charge_minos;
    //scan information
    double vtxx_scan;
    double vtxy_scan;
    double vtxz_scan;
    int neutrino_scan;
    int maybeneutrino_scan;
    int ntracks_scan;
    int nshowers_scan; 
    
    //mctruth information
    int nuPDG_truth;
    int ccnc_truth;
    int mode_truth;
    double enu_truth;
    double Q2_truth;
    double W_truth;
    double nuvtxx_truth;
    double nuvtxy_truth;
    double nuvtxz_truth;
    double lep_mom_truth;
    double lep_dcosx_truth;
    double lep_dcosy_truth;
    double lep_dcosz_truth;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel; 
    std::string fGenieGenModuleLabel;
    std::string fKingaModuleLabel;
    std::string fLineMergerModuleLabel;
    std::string fDbscanModuleLabel;
    std::string fClusterModuleLabel;   
    std::string fTrackModuleLabel;
    std::string fEndPoint2DModuleLabel;
    std::string fVertexModuleLabel;
    std::string fMINOSModuleLabel;
    std::string fTrackMatchModuleLabel;
    std::string fScanModuleLabel;
    std::string fPOTModuleLabel;
    double fvertextrackWindow;
    double fvertexclusterWindow;
    double fboundaryWindow;
    
     int no_primaries;
    int *primaries_pdg;
  double *Eng;
  double *Px;
  double *Py;
  double *Pz;
  double *StartPointx;
  double *StartPointy;
  double *StartPointz;
  double *EndPointx;
  double *EndPointy;
  double *EndPointz;
  int *NumberDaughters;
  //from genie:
 int genie_no_primaries;
  double *genie_primaries_pdg;
  double *genie_Eng;
   double *genie_Px;
    double *genie_Py;
     double *genie_Pz;
      double *genie_P;
       int *genie_status_code;
       double *genie_mass;
       int *genie_trackID;
       int *genie_ND;
       int *genie_mother;
       
       int no_hits;
       int *hit_plane;
  int *hit_wire;
  int *hit_channel;
  double *hit_peakT;
  double *hit_charge;
  }; // class CCQEAnalysisTree
}

#endif 
