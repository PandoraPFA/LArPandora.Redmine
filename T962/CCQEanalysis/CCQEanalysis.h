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
    TH1F* proton_track_length;
    TH1F* Vertex_x;
    TH1F* Vertex_y;
    TH1F* Vertex_z;





  }; 

}

#endif 
