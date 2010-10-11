////////////////////////////////////////////////////////////////////////
/// \file  MergeData.cxx
/// \brief Module for argoneut data
///
/// \version $Id: MergeData.cxx,v 1.20 2010/04/30 16:10:37 brebel Exp $
/// \author  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////////

#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "T962_MergeData/MergeData.h"

namespace merge{

  // A macro required for a JobControl module.
  DEFINE_FWK_MODULE(MergeData);

}
