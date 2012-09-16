////////////////////////////////////////////////////////////////////////
//
// MergeData_Paddles class:
// For each ArgoNeuT event, find corresponding Paddles info. and 
// create a Paddles object to store the info. in the event record
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

// LArSoft include
#include "T962/MergeData/MergeDataPaddles.h"

namespace merge{

  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(MergeDataPaddles);

}
