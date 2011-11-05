/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

#include "Geometry/geo.h"

class TH1F;

namespace recob { class Hit; }

namespace cluster{
   
  //--------------------------------------------------------------- 
  class DBScanService {
  public:
    
    
    DBScanService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    virtual ~DBScanService();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void InitScan(art::PtrVector<recob::Hit>& allhits, std::set<unsigned int> badChannels);
    double getSimilarity(const std::vector<double> v1, const std::vector<double> v2); 
    std::vector<unsigned int> findNeighbors( unsigned int pid, double threshold, double threshold2);
    void computeSimilarity();
    void run_cluster();     
    double getSimilarity2(const std::vector<double> v1, const std::vector<double> v2); 
    void computeSimilarity2();
    double getWidthFactor(const std::vector<double> v1, const std::vector<double> v2); 
    void computeWidthFactor();
      

    std::vector<std::vector<unsigned int> > fclusters;               ///< collection of something
    std::vector<std::vector<double> >       fps;                     ///< the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///< mapping point_id -> clusterId     
    std::vector<std::vector<double> >       fsim;                    ///<
    std::vector<std::vector<double> >       fsim2;            	     ///<
    std::vector<std::vector<double> >       fsim3;            	     ///<
     
  private:
      
    // eps radius
    // Two points are neighbors if the distance 
    // between them does not exceed threshold value.
    double fEps;
    double fEps2;
     
    //minimum number of points
    unsigned int fMinPts;
      
    // noise vector
    std::vector<bool>      fnoise;						     
    std::vector<bool>      fvisited;					     
    std::vector<double>    fWirePitch;   ///< the pitch of the wires in each plane
    std::set<unsigned int> fBadChannels; ///< set of bad channels in this detector
  };
} // namespace
