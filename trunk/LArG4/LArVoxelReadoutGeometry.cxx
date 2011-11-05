////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadoutGeometry.cxx
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
///
/// \version $Id: LArVoxelReadoutGeometry.cxx,v 1.3 2009/03/31 17:58:39 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "G4Base/DetectorConstruction.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "LArG4/LArVoxelReadoutGeometry.h"
#include "LArG4/LArVoxelReadout.h"
#include "Simulation/LArVoxelCalculator.h"
#include "Geometry/geo.h"

// G4 includes
#include <G4VUserParallelWorld.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4LogicalVolume.hh>
#include <G4VisAttributes.hh>
#include <G4VSolid.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>
#include <G4VSensitiveDetector.hh>
#include <G4SDManager.hh>
#include <G4Material.hh>
#include <G4String.hh>
#include <globals.hh>

#include <vector>
#include <cmath>

namespace larg4 {

  ////////////////////////////////////////////////////////////////////
  // Declarations needed for the LAr TPC volume search.

  // The lists of the translations and rotations of volumes; built up
  // as we search for the LAr TPC within the world volume.
  typedef std::vector< G4Transform3D > transformList_type;
  transformList_type                   transformList;

  // The LAr TPC volume name, as-is and with an underscore.
  G4String larTPCName;
  G4String larTPCNameUnderscore;

  // The recursive routine that will do the searching.
  G4VPhysicalVolume* volumeSearch( G4VPhysicalVolume* );
  ////////////////////////////////////////////////////////////////////


  // Constructor and destructor.
  LArVoxelReadoutGeometry::LArVoxelReadoutGeometry(const G4String name)
    : G4VUserParallelWorld(name)
  {}

  LArVoxelReadoutGeometry::~LArVoxelReadoutGeometry() 
  {}

  void LArVoxelReadoutGeometry::Construct()

  {
    // With a "parallel geometry", Geant4 has already created a clone
    // of the world physical and logical volumes.  We want to place
    // the LAr TPC, and only the LAr TPC, within this cloned world.

    // Get the parallel world physical volume.
    G4VPhysicalVolume* parallelPhysical = GetWorld();

    // Now we want to place a parallel LAr TPC volume within this
    // parallel world volume.  We only want to duplicate the LAr TPC
    // volume; any other volumes in the "official" geometry are going
    // to be ignored.  Our parallel world will consist only of the
    // world volume and the LAr TPC volume.

    // The hard part is finding the total transformation between the
    // LAr TPC and the world volume, since any number of nested volume
    // may occur between the two in the "official" geometry.

    // The most practical way to find the LAr TPC this is to "hunt"
    // through the volume hierarchy until we find a volume is either
    // the same as or begins with the LAr TPC volume name.

    // Set up the recursive search: We'll need the name of the volume
    // for which we're searching.  Get it from the LArSoft
    // detector-description interface.
    art::ServiceHandle<geo::Geometry> geo;
    larTPCName = static_cast<G4String>( geo->GetLArTPCVolumeName() );

    // Depending on how the GDML file was converted into Geant4, the
    // volume name may either be what we just got, or it may begin
    // with that name following by an underscore ("_").
    larTPCNameUnderscore = larTPCName + "_";

    // Make sure the "chains" of translations and rotations are clear.
    transformList.clear();

    // Look at the original "mass world" of the detector description.
    G4VPhysicalVolume* worldPhysical = g4b::DetectorConstruction::GetWorld();

    // Do the recursive search.
    G4VPhysicalVolume* larTPCPhysical = volumeSearch( worldPhysical );
    if ( larTPCPhysical == 0 ){
      mf::LogError("LArVoxelReadoutGeometry") << "No volume with the name '" 
					      << larTPCName
					      << "' or beginning with '" 
					      << larTPCNameUnderscore
					      << "' was found.  If I can't find the "
					      << "LAr TPC, I can't make voxels.";
      throw cet::exception("LArVoxelReadoutGeometry") << "Can't find LAr TPC volume";
    }
    
    // We have the LAr TPC volume.  We also have the successive
    // geometric transformations (rotations and translations) of all
    // the volumes containing the LAr TPC.  Combine them all into a
    // single transformation within the world volume.

    // Note: A rotation followed by a translation gives different
    // results than a translation followed by a rotation.  I don't
    // know which order Geant4 uses.  However, if I use Geant4's
    // G4Transforms3D (which is really CLHEP::Transform3D) then I
    // don't have to know.

    // Start with a "zero" transform.
    G4Transform3D transform;

    // Go through the list of transforms, in reverse list order
    // (outermost to innermost), and apply them.
    for ( transformList_type::reverse_iterator r = transformList.rbegin();
	  r != transformList.rend(); ++r ){
      transform = transform * (*r);
    }

    // Get the LAr TPC volume, and its shape.
    G4LogicalVolume* larTPCLogical = larTPCPhysical->GetLogicalVolume();
    G4VSolid* larTPCShape = larTPCLogical->GetSolid();

    // We're not going to exactly duplicate the LAr TPC in our
    // parallel world.  We're going to construct a box of voxels.
    // What should the size and position of that box be?

    // To get our first hints, we need the overall dimensions of a
    // "bounding box" that contains the shape.  For now, we'll allow
    // two possible shapes: a box (by the far the most likely) and a
    // cylinder (for bizarre future detectors that I know nothing
    // about).

    G4VSolid* copyTPCShape = 0;
    G4double larTPCHalfXLength = 0;
    G4double larTPCHalfYLength = 0;
    G4double larTPCHalfZLength = 0;
    G4Box* tpcBox = dynamic_cast< G4Box* >( larTPCShape );
    if ( tpcBox != 0 ){
      larTPCHalfXLength = tpcBox->GetXHalfLength();
      larTPCHalfYLength = tpcBox->GetYHalfLength();
      larTPCHalfZLength = tpcBox->GetZHalfLength();
      copyTPCShape = (G4VSolid*) (new G4Box( *tpcBox ));
    }
    else{
      // It's not a box.  Try a cylinder.
      G4Tubs* tube = dynamic_cast< G4Tubs* >( larTPCShape );
      if ( tube != 0 ){
	larTPCHalfXLength = tube->GetOuterRadius();
	larTPCHalfYLength = tube->GetOuterRadius();
	larTPCHalfZLength = tube->GetZHalfLength();
	copyTPCShape = (G4VSolid*) (new G4Tubs( *tube ));
      }
      else{
	mf::LogError("LArVoxelReadoutGeometry") 
	  << "The LAr TPC volume is not a box or a tube. "
	  << "This routine can't convert any other shapes.";
	throw cet::exception("LArVoxelReadoutGeometry") << "Unknown shape in readout geometry";
      }
    }

    LOG_DEBUG("LArVoxelReadoutGeometry")
      << ": larTPCHalfXLength=" << larTPCHalfXLength
      << ": larTPCHalfYLength=" << larTPCHalfYLength
      << ": larTPCHalfZLength=" << larTPCHalfZLength;

    // Get some constants from the LAr voxel information object.
    // Remember, ROOT uses cm.
    art::ServiceHandle<sim::LArVoxelCalculator> lvc;
    G4double voxelSizeX = lvc->VoxelSizeX() * cm;
    G4double voxelSizeY = lvc->VoxelSizeY() * cm;
    G4double voxelSizeZ = lvc->VoxelSizeZ() * cm;
    G4double voxelOffsetX = lvc->VoxelOffsetX() * cm;
    G4double voxelOffsetY = lvc->VoxelOffsetY() * cm;
    G4double voxelOffsetZ = lvc->VoxelOffsetZ() * cm;

    LOG_DEBUG("LArVoxelReadoutGeometry") << ": voxelSizeX=" << voxelSizeX
					 << ", voxelSizeY=" << voxelSizeY
					 << ", voxelSizeZ=" << voxelSizeZ;
    LOG_DEBUG("LArVoxelReadoutGeometry") << ": voxelOffsetX=" << voxelOffsetX
					 << ", voxelOffsetY=" << voxelOffsetY
					 << ", voxelOffsetZ=" << voxelOffsetZ;

    // We want our voxelization region to be an integer multiple of
    // the voxel sizes in all directions; if we didn't do this, we
    // might get into trouble when we start playing with replicas.
    // Compute the the dimensions of our voxelization to be about the
    // size of the LAr TPC region, adjusted to be an integer number of
    // voxels in all directions.

    G4double numberXvoxels = 2.*larTPCHalfXLength / voxelSizeX;
    G4double numberYvoxels = 2.*larTPCHalfYLength / voxelSizeY;
    G4double numberZvoxels = 2.*larTPCHalfZLength / voxelSizeZ;
    numberXvoxels = trunc(numberXvoxels) + 1.;
    numberYvoxels = trunc(numberYvoxels) + 1.;
    numberZvoxels = trunc(numberZvoxels) + 1.;
    G4double voxelBoxHalfX = numberXvoxels * voxelSizeX / 2.;
    G4double voxelBoxHalfY = numberYvoxels * voxelSizeY / 2.;
    G4double voxelBoxHalfZ = numberZvoxels * voxelSizeZ / 2.;

    LOG_DEBUG("LArVoxelReadoutGeometry") << ": voxelBoxHalfX=" << voxelBoxHalfX
					 << ", voxelBoxHalfY=" << voxelBoxHalfY
					 << ", voxelBoxHalfZ=" << voxelBoxHalfZ;

    // Now we have a box that will include an integer number of voxels
    // in each direction.  Note that the material is irrelevant for a
    // "parallel world."
    G4Box* voxelBox = new G4Box("VoxelBox",voxelBoxHalfX,voxelBoxHalfY,voxelBoxHalfZ);
    G4LogicalVolume* voxelBoxLogical = new G4LogicalVolume( voxelBox,
							    0,
							    "VoxelizationLogicalVolume" );

    // If we general an event display within Geant4, we won't want to
    // see this box.
    G4VisAttributes* invisible = new G4VisAttributes();
    invisible->SetVisibility(false);
    voxelBoxLogical->SetVisAttributes(invisible);

    // We have a "box of voxels" that's the right dimensions, but we
    // have to know exactly where to put it.  The user has the option
    // to offset the voxel co-ordinate system.  We want to place our
    // box so the edges of our voxels align with that co-ordinate
    // system.  In effect, we want to offset our "box of voxels" by
    // the user's offsets, modulo the size of the voxel in each
    // direction.

    G4double offsetInVoxelsX = voxelOffsetX / voxelSizeX;
    G4double offsetInVoxelsY = voxelOffsetY / voxelSizeY;
    G4double offsetInVoxelsZ = voxelOffsetZ / voxelSizeZ;
    G4double fractionOffsetX = offsetInVoxelsX - trunc(offsetInVoxelsX);
    G4double fractionOffsetY = offsetInVoxelsY - trunc(offsetInVoxelsY);
    G4double fractionOffsetZ = offsetInVoxelsZ - trunc(offsetInVoxelsZ);
    G4double offsetX = fractionOffsetX * voxelSizeX;
    G4double offsetY = fractionOffsetY * voxelSizeY;
    G4double offsetZ = fractionOffsetZ * voxelSizeZ;

    // Now we know how much to offset the "box of voxels".  Include
    // that in the transformation of the co-ordinates from world
    // volume to LAr TPC volume.
    transform = G4Translate3D( offsetX, offsetY, offsetZ ) * transform;

    LOG_DEBUG("LArVoxelReadoutGeometry") << ": offsetX=" << offsetX
					 << ", offsetY=" << offsetY
					 << ", offsetZ=" << offsetZ;

    LOG_DEBUG("LArVoxelReadoutGeometry") << ": transform = \n";
    for ( G4int i = 0; i < 3; ++i ){
      for ( G4int j = 0; j < 4; ++j ){ 
	LOG_DEBUG("LArVoxelReadoutGeometry") << transform[i][j] << " ";
      }
      LOG_DEBUG("LArVoxelReadoutGeometry") << "\n";
    }

    // Place the box of voxels, with the accumulated transformations
    // computed above.
    new G4PVPlacement( transform,
		       "VoxelizationPhysicalVolume",
		       voxelBoxLogical,
		       parallelPhysical,
		       false,           // Only one volume
		       0);              // Copy number

    // Now we've fill our "box of voxels" with the voxels themselves.
    // We'll do this by sub-dividing the volume in x, then y, then z.

    // Create an "x-slice".
    G4Box* xSlice = new G4Box("xSlice",voxelSizeX/2.,voxelBoxHalfY,voxelBoxHalfZ);
    G4LogicalVolume* xSliceLogical = new G4LogicalVolume( xSlice, 0, "xLArVoxelSlice" );
    xSliceLogical->SetVisAttributes(invisible);

    // Use replication to slice up the "box of voxels" along the x-axis.
    new G4PVReplica( "VoxelSlicesInX",
		     xSliceLogical,
		     voxelBoxLogical,
		     kXAxis,
		     G4int( numberXvoxels ),
		     voxelSizeX );

    // Now do the same thing, dividing that x-slice along the y-axis.
    G4Box* ySlice = new G4Box("ySlice",voxelSizeX/2.,voxelSizeY/2., voxelBoxHalfZ);
    G4LogicalVolume* ySliceLogical = new G4LogicalVolume( ySlice, 0, "yLArVoxelSlice" );
    ySliceLogical->SetVisAttributes(invisible);
    new G4PVReplica( "VoxelSlicesInY",
		     ySliceLogical,
		     xSliceLogical,
		     kYAxis,
		     G4int( numberYvoxels ),
		     voxelSizeY );
    
    // Now divide the y-slice along the z-axis, giving us our actual voxels.
    G4Box* zSlice = new G4Box("zSlice",voxelSizeX/2.,voxelSizeY/2., voxelSizeZ/2.);
    G4LogicalVolume* voxelLogical = new G4LogicalVolume( zSlice, 0, "LArVoxel" );
    voxelLogical->SetVisAttributes(invisible);
    new G4PVReplica( "LArVoxel",
		     voxelLogical,
		     ySliceLogical,
		     kZAxis,
		     G4int( numberZvoxels ),
		     voxelSizeZ );

    // Define the sensitive detector for the voxel readout.  This
    // routine will be called every time a particle deposits energy in
    // a voxel that overlaps the LAr TPC.
    LArVoxelReadout* larVoxelReadout = new LArVoxelReadout();

    // Tell Geant4's sensitive-detector manager that the voxel SD
    // class exists.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    sdManager->AddNewDetector(larVoxelReadout);

    // Set the sensitive detector of the LAr TPC to be the voxel readout.
    voxelLogical->SetSensitiveDetector(larVoxelReadout);
  }

  //---------------------------------------------------------------
  G4VPhysicalVolume* volumeSearch( G4VPhysicalVolume* physicalVolume )
  {
    // Is this volume the LAr TPC?  (The second test means "does
    // volumeName begin with larTPCNameUnderscore?")
    G4String volumeName = physicalVolume->GetName();
    if ( volumeName == larTPCName  ||
	 volumeName.find( larTPCNameUnderscore, 0, larTPCNameUnderscore.length() ) == 0 ){
      // We found it!  Return it to the caller.
      return physicalVolume;
    }

    // This volume isn't the LAr TPC, but maybe one of its daughters
    // is.  Search each daughter.
    G4LogicalVolume* logicalVolume = physicalVolume->GetLogicalVolume();
    G4int numberDaughters = logicalVolume->GetNoDaughters();
    for ( G4int i = 0; i != numberDaughters; ++i ){
      G4VPhysicalVolume* daughter = logicalVolume->GetDaughter(i);

      // Here's the recursive call: See if the daughter contains the
      // LAr TPC.
      G4VPhysicalVolume* result = volumeSearch( daughter );
      
      if ( result != 0 ){
	// We found it!
	
	// Save the translation and rotation for the daughter
	// volume.
	G4ThreeVector translation = daughter->GetObjectTranslation();
	G4RotationMatrix rotation = daughter->GetObjectRotationValue();
	G4Transform3D transform( rotation, translation );
	transformList.push_back( transform );
	
	// Return the successful result to the caller.
	return result;
      }
    } // for each daughter

    // If we get to this point, we've searched all the daughters of
    // the current volume and didn't find the LAr TPC.  Return a
    // failed search.
    return 0;
  }

} // namespace larg4
