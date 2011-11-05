////////////////////////////////////////////////////////////////////////
// $Id: KingaClusterAna.cxx,v 1.36 2010/09/15  kpartyka Exp $
//
// KingaCluster class
//
// 
//
////////////////////////////////////////////////////////////////////////
#ifndef KingaCluster_H
#define KingaCluster_H

#include "TMath.h"
#include <vector>
#include <string>
#include "art/Persistency/Common/PtrVector.h" 


#include "art/Framework/Core/EDProducer.h"

class TH1F;
class TH2F;

namespace recob { class Hit; }
namespace cluster {
   
  class KingaCluster : public art::EDProducer {
    
  public:
    
    explicit KingaCluster(fhicl::ParameterSet const& pset); 
    KingaCluster();
    ~KingaCluster();
         
    void produce(art::Event& evt);
    void beginJob();
    void AngularDistribution(unsigned int tpc, unsigned int plane);
    void FindMax(unsigned int tpc, unsigned int plane);
    //void FinalPeaks();
    void FitAngularDistributions();
    void FindClusters(unsigned int tpc, unsigned int plane);
    void ReassignHitID(unsigned int tpc, unsigned int plane,unsigned int HitPosition,unsigned int WrongPeakNo);
 
  
     
  private:
  
    std::string fDBScanModuleLabel;  
    std::string fGenieGenModuleLabel;  
    std::string fEndPoint2DModuleLabel;
   
    TH1F *fh_theta_ind;
    TH1F *fh_theta_coll;
    TH1F *fh_theta_ind_2D;
    TH1F *fh_theta_coll_2D;
    TH1F *fh_theta_coll_Area;
    TH1F *fh_theta_ind_Area;
    TH1F *Hit_Area_Ind;
    TH1F *Hit_Area_Coll;
   
     art::PtrVector<recob::Hit> allhits;
    std::vector<int> maxBin;    //stores bin # of local maximum
    std::vector<int> MaxStartPoint; //bin no of the starting point of a peak
    std::vector<int> MaxEndPoint;  //bin no of the end point of a peak
    std::vector<int> MaxStartPointTheta; //theta value of the starting point of a peak
    std::vector<int> MaxEndPointTheta; //theta value of the end point of a peak
    std::vector<int> fwire_vertex, fwire_vertex_reco;
    std::vector<double> ftime_vertex, ftime_vertex_reco;
    std::vector<unsigned int> HitsWithClusterID;
    double ftimetick; //get from parameterset
    double fdriftvelocity;  //get from paramtereset 9either k and V)
    double fpi;
    int fMC;
    double MCvertex [3];

    std::vector<double> maxBinValues;
    std::vector<double> OriginalmaxBinValues;
    std::vector<int> SortedMaxBin;
    std::vector<int> FinalPeaks;
    int fpeaks_found; //flag to determine whether the program should continue or not
    bool need_to_reassign_hitsIDs;
    bool go_ahead_at_reassign;
  protected:

    
  };
  

}



#endif // KingaCluster_H
