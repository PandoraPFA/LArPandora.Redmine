////////////////////////////////////////////////////////////////////////
//
// 
// \author kinga.partyka@yale.edu
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef KingaClusterAna_H
#define KingaClusterAna_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


class TH1F;
class TH2F;
namespace recob { class Hit; }
class TTree;
///Cluster finding and building 
namespace cluster {

   
  class KingaClusterAna : public art::EDAnalyzer {

  public:
          
    explicit KingaClusterAna(fhicl::ParameterSet const& pset); 
    virtual ~KingaClusterAna();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

  private:
  double fMCvertex [3];
  std::vector<int> fwire_vertex, fwire_vertex_reco;
  std::vector<double> ftime_vertex, ftime_vertex_reco;
  int fkingaCl_p0,fkingaCl_p1;
  double ftimetick;
  double fdriftvelocity; 
  std::string fKingaModuleLabel;
  std::string fLineMergerModuleLabel;
  std::string fEndPoint2DModuleLabel;
  std::string fClusterCheaterModuleLabel;
  std::string fGenieGenModuleLabel;
  std::string fLArGeantModuleLabel;
  
  
  std::vector<int> fclusters_planeNo_reco_;
  std::vector<double> fStart_pt_w_reco_;
  std::vector<double> fStart_pt_t_reco_;
    int fkingaCl_near_vertex_p0;
    int fkingaCl_near_vertex_p1;
    int fcheatedCl_p0;
    int fcheatedCl_p1;
    int flinemergerCl_p0;
    int flinemergerCl_p1;
    int fcheatedCl_near_vertex_p0;
    int fcheatedCl_near_vertex_p1;
    int flinemergerCl_near_vertex_p0;
    int flinemergerCl_near_vertex_p1;
    
  TH1F *fdiff_time_vtx_p0;
  TH1F *fdiff_wire_vtx_p0;
  TH1F *fdiff_wire_vtx_p1; 
  TH1F *fdiff_no_vertex_clusters_p0;
  TH1F *fdiff_no_vertex_clusters_p1;
  TH1F *fdiff_no_vertex_linemergerclusters_p0;
  TH1F *fdiff_no_vertex_linemergerclusters_p1;
  TH1F *fNoProtonTracks_p0_cheatedCl;
  TH1F *fNoProtonTracks_p1_cheatedCl;
  TH1F *fNoProtonTracks_p0_linemergerCl;
  TH1F *fNoProtonTracks_p1_linemergerCl;
  TH1F *fNoProtonTracks_p0_kingaCl;
  TH1F *fNoProtonTracks_p1_kingaCl;
    
   //TTree:
     TTree* fTree;
     int frun;
    int fevent;
    double ftime_vertex_true;
    int fno_clusters_true;
    int fno_clusters_reco;
    int fno_clusters_linemerger;
   
    double *fwire_vertex_true;
    double *fTTree_wire_vertex_reco;
    double *fTTree_time_vertex_reco;
    int *fclusters_planeNo_true;
    int *fclusters_planeNo_reco;
    double *fStart_pt_w_true;
    double *fStart_pt_t_true;
    double *fStart_pt_w_reco;
    double *fStart_pt_t_reco;
    double *fStart_pt_t_linemerger;
    double *fStart_pt_w_linemerger;
     int *flinemergerclusters_planeNo;
    int *fcheated_cluster_size;
    int *flinemerger_cluster_size;




  }; // class KingaClusterAna

}

#endif 
