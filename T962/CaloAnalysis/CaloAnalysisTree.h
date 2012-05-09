////////////////////////////////////////////////////////////////////////
// 
//
// kinga.partyka@yale.edu
// 
////////////////////////////////////////////////////////////////////////
#ifndef CaloAnalysisTree_H
#define CaloAnalysisTree_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "RecoBase/Track.h"

#include <vector>
#include <string>

class TH1F;
class TH2D;
class TTree;

///Track finding and building 
namespace t962 {

   
  class CaloAnalysisTree : public art::EDAnalyzer {

  public:
          
    explicit CaloAnalysisTree(fhicl::ParameterSet const& pset); 
    virtual ~CaloAnalysisTree();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();
    void beginSubRun(const art::SubRun& sr);

    
    
    
  protected:
  
 
  
  private:
    double LifetimeCorrection(float time);
    void ResetVars();
    double fMCvertex [3];
    
    TTree* fTree;
    //run information
    TH1F *diff_length_reco_truth;
    TH1F *diff_KE_reco_truth;
    TH1F *diff_length_truth1_vs_truth2;
    TH1F *diff_length_reco_truth_1geant;
    TH1F *diff_KE_reco_truth_1geant;
    TH1F *diff_length_truth1_vs_truth2_1geant;
    
    double *trk_length_truth;
    double *trk_length_straight_line_truth;
    double *trk_length_reco;
    double *Kin_Eng_reco;
    double Kin_Eng_truth;
    double KE_sum_;
    int run;
    int event;
    double pot;
    int isdata;
    //reconstructed information    
    double vtxx_reco;
    double vtxy_reco;
    double vtxz_reco;
   
    int ntracks_reco;         //number of reconstructed tracks
   
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
    
    //mctruth information
   
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
    
     int no_geant_particles;
     double *twodvtx_w_reco;
    double *twodvtx_t_reco;
    double *twodvtx_w_truth;
    double *twodvtx_t_truth;
    
    int *pdg;
    int *process_primary;
  double *Eng;
  double *Eng_at_endtrack;
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
 
 
  double fTrkPitchC;
  int fnhitsCOL; 
  
  
  }; // class CaloAnalysisTree
}

#endif 
