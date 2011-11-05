////////////////////////////////////////////////////////////////////////
//
// LArTracker class
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef LARTRACKER_H
#define LARTRACKER_H

#include "art/Framework/Core/EDProducer.h" 

#include <vector>
#include <string>

#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h"
#include "RecoBase/recobase.h"

///track finding
namespace trkf {
   
   class LArTracker : public art::EDProducer {
    
   public:
    
      explicit LArTracker(fhicl::ParameterSet const& pset);
      ~LArTracker();
    
      //////////////////////////////////////////////////////////
      void reconfigure(fhicl::ParameterSet const& p);
      void beginJob();
      void endJob();
      void produce(art::Event& evt); 

   private:

      //Return the position along the drift ("x") direction, 
      //adjusting for spacing between planes. 
      double DriftCoordinate(int plane, double time);

      //Do two input clusters have endpoints that are consistent in time?
      bool   ClusterEndPointsMatch(art::Ptr<recob::Cluster> c1, art::Ptr<recob::Cluster> c2); 

      //Helper functions to avoid double-counting
      bool   MultiPlane  (art::PtrVector<recob::Hit> hitList);
      bool   MultiPlane  (art::PtrVector<recob::Cluster> clusterList);
      bool   FindCluster (art::PtrVector<recob::Cluster> clusterList, art::Ptr<recob::Cluster> c);
      bool   FindHit     (art::PtrVector<recob::Hit> hitList,         art::Ptr<recob::Hit> h);
      bool   FindHit     (std::vector<recob::SpacePoint> spList,      art::Ptr<recob::Hit> h);
      bool   FindHit     (recob::Track track,                         art::Ptr<recob::Hit> h);
         
      std::string     fClusterModuleLabel;// label for input cluster collection
      double          ftmatch; // tolerance for time matching (in cm)
   protected: 
    
  
   }; // class LArTracker

}

#endif // LARTRACKER_H
