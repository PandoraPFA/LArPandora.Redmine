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
#include "T962/T962_Objects/Paddles.h"
#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"
#include "T962/T962_Objects/ScanInfo.h"
//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//
template class std::vector<t962::MINOS>;
template class std::vector<t962::MINOSTrackMatch>;
template class std::vector<t962::Paddles>;
template class std::vector<t962::ScanInfo>;

template class art::Wrapper< std::vector<t962::MINOS>     >;
template class art::Wrapper< std::vector<t962::MINOSTrackMatch> >;
template class art::Wrapper< std::vector<t962::Paddles>   >;
template class art::Wrapper< std::vector<t962::ScanInfo>   >;
