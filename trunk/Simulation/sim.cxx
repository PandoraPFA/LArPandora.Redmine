////////////////////////////////////////////////////////////////////////
/// \file  sim.cxx
/// \brief useful tools for simulation
///
/// 
///
/// \version $Id: evgenbase.cxx,v 1.1 2011/06/29 17:13:40 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Simulation/sim.h"

#include "TRandom3.h"

namespace sim{

  unsigned int GetRandomNumberSeed()
  {

    // the maximum allowed seed for the art::RandomNumberGenerator
    // is 900000000. Use TRandom3 to get the seed value in that range.
    // Instantiating TRandom3 with a 0 means that its seed is set based
    // on the TUUID and should always be random, even for jobs running on the
    // same machine
    TRandom3 rand(0);
    return rand.Integer(900000000);
  }

}
