////////////////////////////////////////////////////////////////////////
// $Id: SimWireT962.cxx,v 1.22 2010/04/23 20:30:53 seligman Exp $
//
// SimWireT962 class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////
// LArSoft includes
#include "DetSim/SimWireT962.h"

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace detsim{

  DEFINE_ART_MODULE(SimWireT962);

}
