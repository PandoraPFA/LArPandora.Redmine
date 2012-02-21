////////////////////////////////////////////////////////////////////////
//
// MergeData_Beam class:
// For each ArgoNeuT event, find corresponding NuMI spill info. and 
// create a BeamInfo object to store the info. in the event record
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

// LArSoft include
#include "T962/MergeData/MergeDataBeam.h"

namespace merge{

  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(MergeDataBeam);

}
