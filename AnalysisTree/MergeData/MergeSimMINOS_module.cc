////////////////////////////////////////////////////////////////////////
//
// MergeDataMINOS class:
// For each ArgoNeuT event, find corresponding MINOS info. using timing only, 
// and create a MINOS object to store the info. in the event record
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

// LArSoft includes
#include "T962/MergeData/MergeSimMINOS.h"

namespace merge{

  DEFINE_ART_MODULE(MergeSimMINOS);

}
