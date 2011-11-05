////////////////////////////////////////////////////////////////////////////
// \version $Id: Shower.cxx,v 1.2 2010/02/15 20:32:46 brebel Exp $
//
// \brief Definition of shower object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#include "RecoBase/Shower.h"
#include <string>
#include <iostream>

namespace recob{

  //----------------------------------------------------------------------
  Shower::Shower() :
    Prong()
  {
  }

  //----------------------------------------------------------------------
  Shower::Shower(art::PtrVector<recob::Cluster> &clusters,  std::vector<recob::SpacePoint> spacepoints) :
    Prong(clusters,spacepoints)
  {
  }

  //----------------------------------------------------------------------
  Shower::~Shower()
  {
  }

}
