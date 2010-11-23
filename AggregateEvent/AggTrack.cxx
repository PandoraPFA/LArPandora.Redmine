////////////////////////////////////////////////////////////////////////
// $Id: AggTrack.cxx,v 1.1 2010/09/02 17:25:11 echurch Exp $
//
// AggTrack class
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

#include "AggregateEvent/AggTrack.h"
#include "TH1.h"

namespace aggr {

  AggTrack::AggTrack(edm::Ptr<recob::Track> dum) :
    copiedTrack(dum)
  {
  }

  edm::Ptr<recob::Track> AggTrack::GetTrack()
  {
    edm::Ptr<recob::Track> t;
    t = copiedTrack;
    return t;
  }

  edm::PtrVector<recob::Track> AggTrack::GetShowers()
  {
    edm::PtrVector<recob::Shower> slist;
    for (int ii = 0 ; ii<matchedShowers.size(); ii++)
      {
	slist.push_back(matchedShowers[ii]);
      }
    return slist;
  }

}
