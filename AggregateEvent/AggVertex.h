
#ifndef AGGVTX_H
#define AGGVTX_H

#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include <vector>
#include <string>

namespace aggr {


  class AggVertex 
  {

  public:

    AggVertex();
    AggVertex(edm::Ptr<recob::Vertex> );
    ~AggVertex();

    void SetMatchedTracks(edm::PtrVector<recob::Track>  dumMatchedTracks){matchedTracks=dumMatchedTracks;};
    edm::Ptr<recob::Vertex> GetVert();
    edm::PtrVector<recob::Track> GetTracks();

  private:

    edm::Ptr<recob::Vertex> copiedVert;
    edm::PtrVector<recob::Track> matchedTracks;
    friend class AggregateVertex;
    friend class AggregateVertexAna;

  }; // class AggVertex

}  // Namespace aggr

#endif // AGGREGATE_H
