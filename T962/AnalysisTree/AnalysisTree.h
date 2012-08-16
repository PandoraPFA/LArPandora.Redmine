////////////////////////////////////////////////////////////////////////
// \version $Id$
// 
// \author tjyang@fnal.gov
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "RecoBase/recobase.h"

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TTree;

const int kMaxTrack = 1000; //maximum number of tracks

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
    
    void HitsPurity(std::vector< art::Ptr<recob::Hit> > const& hits, int& trackid, double& purity);
    void ResetVars();

    TTree* fTree;
    //run information
    int run;
    int event;
    double evttime;
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
    double trkrange[kMaxTrack];
    int    trkpid[kMaxTrack];
    int    trkpidndf[kMaxTrack];
    double trkpidchi2[kMaxTrack];
    double trkmissinge[kMaxTrack];
    double trkmissingeavg[kMaxTrack];
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
    //scan information
    double vtxx_scan;
    double vtxy_scan;
    double vtxz_scan;
    int neutrino_scan;
    int maybeneutrino_scan;
    int ntracks_scan;
    int nshowers_scan; 
    
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
    double lep_mom_truth;
    double lep_dcosx_truth;
    double lep_dcosy_truth;
    double lep_dcosz_truth;
    double beamwgt;
    
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
    int trkcontained_minos;
    int test_charge_minos;
    int hitnuc_truth;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel; 
    std::string fGenieGenModuleLabel;
    std::string fG4ModuleLabel;
    std::string fClusterModuleLabel;   
    std::string fTrackModuleLabel;
    std::string fEndPoint2DModuleLabel;
    std::string fVertexModuleLabel;
    std::string fMINOSModuleLabel;
    std::string fTrackMatchModuleLabel;
    std::string fScanModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    double fvertextrackWindow;
    double fvertexclusterWindow;
    double fboundaryWindow;

    TH1D *hBeamWeight_numu_numode;

  }; // class AnalysisTree
}

#endif 
