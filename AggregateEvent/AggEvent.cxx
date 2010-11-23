////////////////////////////////////////////////////////////////////////
// $Id: AggEvent.cxx,v 1.1 2010/09/02 17:25:11 echurch Exp $
//
// AggEvent class
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

#include "AggregateEvent/AggEvent.h"
#include "TH1.h"

//-----------------------------------------------
aggr::AggEvent::AggEvent(edm::PtrVector<aggr::AggVertex> &av,edm::PtrVector<aggr::AggTrack> &at)
{
}

aggr::AggEvent::~AggEvent() 
{
}


;
