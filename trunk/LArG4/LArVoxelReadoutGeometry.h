////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadoutGeometry.h
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
///
/// \version $Id: LArVoxelReadoutGeometry.h,v 1.1.1.1 2009/02/23 17:28:35 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// This class defines the parallel geometry that will be divided into
/// the three-dimensional voxels for the detector read-out.
/// 
/// Why define a parallel geometry?  Here are some reasons:
///
/// - The regular LAr TPC is one large volume of liquid argon.  When
///   Geant4 does its physics modeling, it can be unconstrained in
///   step size by the voxels.  Only for readout would the steps be
///   sub-divided.
///
/// - There may be more than one kind of readout, depending on a
///   detector's instrumentation (e.g., PMTs in addition to the wire
///   planes).  It's possible that the voxelization appropriate for
///   the wire planes may not be an appropriate readout for the other
///   readouts.  Geant4 allows the construction of multiple parallel
///   readouts, so this mechanism is relatively easy to extend for
///   each type of readout.

#ifndef LArG4_LArVoxelReadoutGeometry_h
#define LArG4_LArVoxelReadoutGeometry_h

#include <G4VUserParallelWorld.hh>
#include <G4String.hh>

namespace larg4 {

  class LArVoxelReadoutGeometry : public G4VUserParallelWorld
  {
  public:
    /// Constructor and destructor.
    LArVoxelReadoutGeometry( const G4String name = "LArVoxelReadoutGeometry" );
    virtual ~LArVoxelReadoutGeometry();

    /// The key method in this class; creates a parallel world view of
    /// those volumes relevant to the LAr voxel readout.  Required of
    /// any class that inherits from G4VUserParallelWorld
    virtual void Construct();
  };

} // namespace larg4

#endif // LArG4_LArVoxelReadoutGeometry_h
