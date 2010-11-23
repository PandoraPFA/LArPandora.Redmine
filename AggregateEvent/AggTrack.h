
#ifndef AGGTRK_H
#define AGGTRK_H

#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include <vector>
#include <string>

namespace aggr {


  class AggTrack 
  {

  public:

    AggTrack();
    AggTrack(edm::Ptr<recob::Track> );
    ~AggTrack();

    void SetMatchedShowers(std::vector<const recob::Track*> );
    edm::Ptr<recob::Track> GetTrack();
    edm::PtrVector<recob::Track> GetShowers();

  private:

    edm::Ptr<recob::Track> copiedTrack;
    edm::PtrVector<recob::Shower> matchedShower;

  }; // class AggTrack

}  // Namespace aggr

#endif // AGGREGATE_H
