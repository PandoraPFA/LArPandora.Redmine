////////////////////////////////////////////////////////////////////////
// $Id: LArFFT.cxx,v 1.11 2010/07/02 20:33:09 bpage Exp $
//
// \file LArFFT_plugin
//
//  This class simplifies implementation of Fourier transforms. 
//  Because all data inputs and outputs are purely real,  the 
//  transforms implemented in this way get a substantial performance
//  increase ~2x.
//
// \author pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "Utilities/LArFFT.h"

namespace util{
 
  DEFINE_ART_SERVICE(LArFFT);

}
