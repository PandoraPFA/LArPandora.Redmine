////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parmeters.cxx
/// \brief Store parameters for running LArG4
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This class exists to pass parameters loaded by the job controlling object to other
// related classes. To load or set a parameter, create an instance of LArG4Parameters and
// use the getting and setting methods,all of which point to one particular instance of
// the class, TheLArG4Parameters
//
// Note - I plan to come back to this class and make it more general (ie load the whole
// config and not use specific get / set methods) soon.
//
// Ben Jones, MIT, March 2010


// Framework includes
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "Simulation/LArG4Parameters.h"

namespace sim {

  DEFINE_ART_SERVICE(LArG4Parameters);

}
