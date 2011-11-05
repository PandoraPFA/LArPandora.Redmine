////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelCalculator_plugin.cc
/// \brief Encapsulates calculation of LArVoxelID and LArVoxel parameters
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

/// This service encapsulates the calculations associated with
/// computing the LArVoxelID, and provides access to the any LArVoxel
/// parameters from the input file(s).

/// Definition: "Voxels" are three-dimensional "pixels"; basically
/// they divide the energy deposition in the LAr into (x,y,z) cubes.
/// Well, hyper-cubes actually, since we have to potentially include
/// divisions in time as well.

// Framework includes
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "Simulation/LArVoxelCalculator.h"

namespace sim {

  DEFINE_ART_SERVICE(LArVoxelCalculator);


} // namespace util
 
