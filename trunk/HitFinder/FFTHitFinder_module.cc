////////////////////////////////////////////////////////////////////////
//
// FFTHitFinder class
//
// pagebri3@msu.edu
//
//  This algorithm is designed to find hits on wires after deconvolution
//  with an average shape used as the input response.
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 

#include "HitFinder/FFTHitFinder.h"

namespace hit{

  DEFINE_ART_MODULE(FFTHitFinder);

} // end of hit namespace
