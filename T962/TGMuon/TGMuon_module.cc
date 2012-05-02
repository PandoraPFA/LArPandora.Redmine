////////////////////////////////////////////////////////////////////////
/// \file  SingleGen_plugin.cc
/// \brief Generator for cosmic-rays
///
/// Module designed to produce a set list of particles for a MC event
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

// lar includes
#include "T962/TGMuon/TGMuon.h"

namespace t962{

  DEFINE_ART_MODULE(TGMuon);

}//end namespace evgen
////////////////////////////////////////////////////////////////////////
