////////////////////////////////////////////////////////////////////////
//
// 
// \author saima@ksu.edu
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef CCRECOEFFTREE_H
#define CCRECOEFFTREE_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "RecoBase/Track.h"

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TTree;

///Track finding and building 
namespace t962 {

   
  class CCRecoEffTree : public art::EDAnalyzer {

  public:
          
    explicit CCRecoEffTree(fhicl::ParameterSet const& pset); 
    virtual ~CCRecoEffTree();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);

    
    
    
  protected:
  
  bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);
  
  private:
    
    void ResetVars();

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
    int hitnuc_truth;    
    int test_charge_minos;

    int nminos_tracks;
    double *trk_charge_minos_all;

    int muon_reco;
    int minos_enter_true;

    double trackstart_x_reco_muon;
    double trackstart_y_reco_muon;
    double trackstart_z_reco_muon;
    double trackexit_x_reco_muon;
    double trackexit_y_reco_muon;
    double trackexit_z_reco_muon; 

    int muon_exits;

    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fVertexModuleLabel;
    std::string fMINOSModuleLabel;
    std::string fTrackMatchModuleLabel;
    std::string fPOTModuleLabel;
    double fboundaryWindow;
  }; // class CCRecoEffTree
}

#endif 
