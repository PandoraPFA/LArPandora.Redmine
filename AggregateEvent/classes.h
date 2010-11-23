//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  echurch@fnal.gov$
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke,
//
// Notes:
// 1) The system is not able to deal with
//    edm::Wrapper<std::vector<std::string> >;
//    The problem is somewhere inside root's reflex mechanism
//    and Philippe Canal says that it is ( as of March 2010) a
//    known problem.  He also says that they do not have any
//    plans to fix it soon.  We can always work around it 
//    by putting the string inside another object.

#include "DataFormats/Common/interface/Wrapper.h"

// nutools includes
#include "AggregateEvent/AggVertex.h"
#include "AggregateEvent/AggTrack.h"
#include "AggregateEvent/AggEvent.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class edm::Wrapper< std::vector<aggr::AggVertex>  >;
template class edm::Wrapper< std::vector<aggr::AggTrack>   >;
template class edm::Wrapper< std::vector<aggr::AggEvent>   >;
