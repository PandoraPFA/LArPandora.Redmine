////////////////////////////////////////////////////////////////////////
//
// 
// \author kinga.partyka@yale.edu
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef DBCLUSTERANA_H
#define DBCLUSTERANA_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


class TH1F;
class TH2F;
///Cluster finding and building 
namespace cluster {

   
  class DBclusterAna : public art::EDAnalyzer {

  public:
          
    explicit DBclusterAna(fhicl::ParameterSet const& pset); 
    virtual ~DBclusterAna();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

  private:
    TH1F* fNoParticles_pdg;
    TH1F* fNoParticles_trackid; 
    TH1F* fNoParticles_trackid_mother;
    TH1F* fNoParticles_trackid_per_event;  
    TH1F* fNoParticles_pdg_per_event;
    TH1F* fCl_for_Muon;
   /*  TH1F* fCl_for_Electron; */
   /*  TH1F* fCl_for_Positron; */
   /*  TH1F* fCl_for_Pion_111; */
   /*  TH1F* fCl_for_Pion_211; */
   /*  TH1F* fCl_for_Pion_m211;*/
   /*  TH1F* fCl_for_Proton;   */
    TH1F* fNoClustersInEvent;
    TH1F* fPercentNoise;
    TH1F* fno_of_clusters_per_track;
    TH1F* fPercent_lost_muon_hits;
    TH1F* fPercent_lost_electron_hits;
    TH1F* fPercent_lost_positron_hits;
    TH1F* fPercent_lost_111_hits;
    TH1F* fPercent_lost_211_hits;
    TH1F* fPercent_lost_m211_hits;
    TH1F* fPercent_lost_2212_hits;
    TH1F* fPercent_lost_2112_hits;

    TH1F* fPercent_lost_muon_energy;
    TH1F* fPercent_lost_electron_energy;
    TH1F* fPercent_lost_positron_energy;
    TH1F* fPercent_lost_111_energy;
    TH1F* fPercent_lost_211_energy;
    TH1F* fPercent_lost_m211_energy;
    TH1F* fPercent_lost_2212_energy;
    TH1F* fPercent_lost_2112_energy;
    TH1F* fEnergy;
    TH2F* fbrian_in;
    TH2F* fbrian_coll;
    
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fClusterFinderModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;

    
      	 
  }; // class DBclusterAna

}

#endif 
