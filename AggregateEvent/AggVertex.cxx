////////////////////////////////////////////////////////////////////////
// $Id: AggVertex.cxx,v 1.1 2010/09/02 17:25:11 echurch Exp $
//
// AggVertex class
//
//  This class will contain the already-discovered vertices
//  and then a vector of track pointers  which we associate
//  to them. The association is done in the AggregateEvent package.
//
// echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "AggregateEvent/AggVertex.h"
#include "TH1.h"

namespace aggr {

  AggVertex::AggVertex(const recob::Vertex* dum) :
    copiedVert(dum)
  {
  }

  edm::Ptr<recob::Vertex> AggVertex::GetVertex()
  {
    edm::Ptr<recob::Vertex> v;
    v = copiedVert;
    return v;
  }

  edm::PtrVector<recob::Track> AggVertex::GetTracks()
  {
    edm::PtrVector<recob::Track> tlist;
    for (int ii = 0 ; ii<matchedTracks.size(); ii++)
      {
	tlist.push_back(matchedTracks[ii]);
      }
    return tlist;
  }

}
