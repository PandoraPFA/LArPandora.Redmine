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

#include "RecoBase/recobase.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class std::vector<recob::Cluster>;
template class std::vector<recob::Hit>;
template class std::vector<recob::Prong>;
template class std::vector<recob::Shower>;
template class std::vector<recob::Track>;
template class std::vector<recob::EndPoint2D>;
template class std::vector<recob::Wire>;
template class std::vector<recob::SpacePoint>;
template class std::vector<recob::Vertex>;
template class std::vector<recob::Event>;

template class art::Wrapper< std::vector<recob::Cluster>    >;
template class art::Wrapper< std::vector<recob::Hit>        >;
template class art::Wrapper< std::vector<recob::Prong>      >;
template class art::Wrapper< std::vector<recob::Shower>     >;
template class art::Wrapper< std::vector<recob::Track>      >;
template class art::Wrapper< std::vector<recob::EndPoint2D> >;
template class art::Wrapper< std::vector<recob::Wire>       >;
template class art::Wrapper< std::vector<recob::SpacePoint> >;
template class art::Wrapper< std::vector<recob::Vertex>     >;
template class art::Wrapper< std::vector<recob::Event>      >;
