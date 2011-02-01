//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by klg
//
// Notes:
// 1) The system is not able to deal with
//    art::Wrapper<std::vector<std::string> >;
//    The problem is somewhere inside root's reflex mechanism
//    and Philippe Canal says that it is ( as of March 2010) a
//    known problem.  He also says that they do not have any
//    plans to fix it soon.  We can always work around it 
//    by putting the string inside another object.

#include "art/Persistency/Common/Wrapper.h"

// nutools includes
#include "T962_MergeData/ScanInfo.h"
#include "T962_MergeData/Paddles.h"
#include "T962_MergeData/MINOS.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class std::vector< std::vector<int> >;
template class std::vector<merge::ScanInfo>;
template class std::vector<merge::MINOS>;
template class std::vector<merge::Paddles>;

template class art::Wrapper< std::vector< std::vector<int> > >;
template class art::Wrapper< std::vector<merge::ScanInfo>  >;
template class art::Wrapper< std::vector<merge::MINOS>     >;
template class art::Wrapper< std::vector<merge::Paddles>   >;
