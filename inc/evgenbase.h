#ifndef EVGENBASE_H
#define EVGENBASE_H
///

#include "EventGeneratorBase/inc/GENIEHelper.h"
#include "EventGeneratorBase/inc/CRYHelper.h"

/// Physics generators for neutrinos, cosmic rays, and others
namespace evgb {
  /// Enumerate mother codes for primary particles. 
  ///
  /// Normally the mother code for a primary particle would be set to
  /// some arbitrary invalid value like -1, however, we can use this
  /// to mark the source of the particle as being either, eg.,
  /// neutrino induced or from cosmic-rays.
  enum {
    kNeutrinoGenerator  = -100,
    kCosmicRayGenerator = -200
  };
}
#endif
