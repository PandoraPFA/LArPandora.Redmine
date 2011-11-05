////////////////////////////////////////////////////////////////////////
/// \file PMTReadoutGeometry.cxx
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
#include "G4Base/DetectorConstruction.h"

#include "LArG4/PMTReadoutGeometry.h"
#include "LArG4/PMTLookup.h"
#include "LArG4/PMTSensitiveDetector.h"
#include "G4PVPlacement.hh"
#include "G4VSolid.hh"
#include "G4SDManager.hh"
#include "sstream"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {

  PMTReadoutGeometry::PMTReadoutGeometry(G4String PMTSensitiveName, const G4String name) :
    G4VUserParallelWorld(name)
  {
    fPMTSensitiveName = PMTSensitiveName;   
  }

  PMTReadoutGeometry::~PMTReadoutGeometry()
  {}

  void PMTReadoutGeometry::Construct()
  {
    mf::LogInfo("PMTReadoutGeometry") << "constructing parallel world, looking for " 
				      << fPMTSensitiveName;

    // Get an empty parallel world volume
    G4VPhysicalVolume * ParallelWorld = GetWorld();
    
    // Start with empty vectors
    PMTVolumes.clear();
    PMTTransformations.clear();

    // Get the primary world volume
    G4VPhysicalVolume * WorldPhysical = g4b::DetectorConstruction::GetWorld();

    // Find the PMT volumes
    std::vector<G4Transform3D> EmptyVector;
    EmptyVector.clear();
    FindVolumes(WorldPhysical, fPMTSensitiveName, EmptyVector );

    // Get the PMT Lookup Table
    PMTLookup * ThePMTLookup = PMTLookup::Instance();

    // Create sensitive detector
    PMTSensitiveDetector * TheSD = new PMTSensitiveDetector("PMTSensitiveDetector");


    if(PMTVolumes.size()>0)
      {

	
	
	// Make placements
	for(unsigned int i=0; i!=PMTVolumes.size(); i++)
	{
	  std::stringstream VolumeName;
	  VolumeName.flush();
	  VolumeName.str("PMTVolume_");
	  VolumeName<<i;
	  
	  G4Transform3D      TheTransform = PMTTransformations.at(i);

	  G4VSolid * TheSolid = PMTVolumes.at(i)->GetSolid();
	  G4Material * TheMaterial = PMTVolumes.at(i)->GetMaterial();
	  G4LogicalVolume * TheLogVolume = new G4LogicalVolume(TheSolid,TheMaterial,"VolumeName.str().c_str()");
	  
	  TheLogVolume->SetSensitiveDetector(TheSD   );
	  

	  G4PVPlacement * ThePlacement 
	    = new G4PVPlacement( TheTransform,
				 VolumeName.str().c_str(),
				 TheLogVolume,
				 ParallelWorld,
				 false,          
				 0);             
	  
	  CLHEP::Hep3Vector trans = ThePlacement->GetTranslation();

	  

	  ThePMTLookup->AddPhysicalVolume(ThePlacement);
	 

	}
      }



  }


  void PMTReadoutGeometry::FindVolumes(G4VPhysicalVolume * PhysicalVolume, G4String PMTName, std::vector<G4Transform3D> TransformSoFar)
  {
    
    // Add the next layer of transformation to the vector
    G4ThreeVector Translation = PhysicalVolume->GetObjectTranslation();
    G4RotationMatrix Rotation = PhysicalVolume->GetObjectRotationValue();
    G4Transform3D NextTransform( Rotation, Translation );
    
    TransformSoFar.push_back(NextTransform);


    // Check if this volume is a PMT
    G4String PMTNameUnderscore = PMTName+"_";
    G4String VolumeName = PhysicalVolume->GetName();   
    if( ( VolumeName == PMTName ) ||
	( VolumeName.find( PMTNameUnderscore,0,PMTNameUnderscore.length() )==0 )
	)
      {

	// We found a PMT! Store its volume and global transformation
	G4ThreeVector     Trans(0,0,0);
	G4RotationMatrix  Rot(0,0,0);
	G4Transform3D TotalTransform(Rot,Trans);
	//for ( std::vector<G4Transform3D>::reverse_iterator it = TransformSoFar.rbegin();
	//     it != TransformSoFar.rend(); ++it )
	


	for ( std::vector<G4Transform3D>::iterator it = TransformSoFar.begin();
	      it != TransformSoFar.end(); ++it )
	  {
	    CLHEP::Hep3Vector trans = (*it).getTranslation();
	    CLHEP::HepRotation rot =  (*it).getRotation();
	    TotalTransform =  TotalTransform * (*it);
	  }

	PMTVolumes.push_back(PhysicalVolume->GetLogicalVolume());	
	PMTTransformations.push_back(TotalTransform);

      }
    else
      {
	// We did not find a PMT.  Keep looking through daughters.
	G4LogicalVolume * LogicalVolume = PhysicalVolume->GetLogicalVolume();
		
	// Loop through the daughters of the volume
	G4int NumberDaughters = LogicalVolume->GetNoDaughters();
	for(G4int i=0; i!=NumberDaughters; ++i)
	  {
	    // Get the ith daughter volume
	    G4VPhysicalVolume * Daughter = LogicalVolume->GetDaughter(i);
	    
	    // Recursively step into this volume
	    FindVolumes(Daughter, PMTName, TransformSoFar);
	  }
      }
  }
    
}
