////////////////////////////////////////////////////////////////////////
//
// 
// \author kinga.partyka@yale.edu
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef CCQEanalysis_H
#define CCQEanalysis_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


class TH1F;
class TH2F;
class TTree;
namespace recob { class Hit; }
///Cluster finding and building 
namespace t962 {

   
  class CCQEanalysis : public art::EDAnalyzer {

  public:
          
    explicit CCQEanalysis(fhicl::ParameterSet const& pset); 
    virtual ~CCQEanalysis();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

  private:
  
  
  //void ResetVars();
  
    TTree* fTree;
    //run information
    int frun;
    int fevent;
    int fbeam;
    double pot;
    int isdata;
     //mctruth information
    int ccnc_truth;
    int mode_truth;
    double Mu_theta_truth;
    double Mu_phi_truth;
    int Mu_charge_truth;
     double lep_dcosx_truth;
    double lep_dcosy_truth;
    double lep_dcosz_truth;
    double lep_mom_truth;
    //reconstructed information 
    
    int fno_primaries;
  //std::vector<int> fprimaries_pdg;
  //std::vector<double> fEng;
  
  int *fprimaries_pdg;
  double *fEng;
  double *fPx;
  double *fPy;
  double *fPz;
  double *fStartPointx;
  double *fStartPointy;
  double *fStartPointz;
  double *fEndPointx;
  double *fEndPointy;
  double *fEndPointz;
//     
  std::string fGenieGenModuleLabel;
  std::string fLArG4ModuleLabel;
  std::string fHitsModuleLabel;
  std::string fClusterFinderModuleLabel;
   
    
    
    TH1F* Mu_theta;
    TH1F* Mu_phi;
    TH1F* mu_charge;
    TH1F* No_protons_in_event;
    TH1F* No_particles_in_event;
    TH1F* P_containment;
    TH1F* P_containment_no_cut;
    TH1F* proton_track_length;
    TH1F* track_l_vs_containment;
    TH1F* p_mom_vs_containment;
    TH1F* track_l_vs_mom;
    
    TH1F* Vertex_x;
    TH1F* Vertex_y;
    TH1F* Vertex_z;





  }; 

}

#endif 
