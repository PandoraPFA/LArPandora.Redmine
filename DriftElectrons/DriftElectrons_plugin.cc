////////////////////////////////////////////////////////////////////////
/// \file  DriftElectrons.cxx
/// \brief Module to drift ionization electrons to wire planes
///
/// \version $Id: DriftElectrons.cxx,v 1.20 2010/04/30 16:10:37 brebel Exp $
/// \author  baller@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DriftElectrons/DriftElectrons.h"

namespace dfe{

  // A macro required for a JobControl module.
  DEFINE_FWK_MODULE(DriftElectrons);

}
