////////////////////////////////////////////////////////////////////////
//
// ScanFilter class:
// Tells the downstream rreconstruction to not process events that
// do not pass certain hand scan criteria
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
 
// LArSoft includes
#include "T962/ScanFilter/ScanFilter.h"

namespace filt {

  DEFINE_ART_MODULE(ScanFilter);

} //namespace filt
