////////////////////////////////////////////////////////////////////////
/// \file PMTReadoutGeometry.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// PMTReadoutGeometry
//
// Ben Jones, MIT, 07/16/10
//
// PMT's are defined with a particular geometry, and only some part of the PMT
// is sensitive.  When several are placed, only the mother volume, vol_PMT is 
// replicated.  Each constituent volume, including the sensitive volume, only
// exists once in the volume store, as a daughter of the mother volume.
//
// Hence to know which PMT a photon steps into when it is inside a sensitive
// volume, we must define a readout geometry to identify which volumes are
// contained within each PMT.  
//
// This class is heavily based on LArG4ReadoutGeometry by Bill Seligman,
// which is very well commented.  See that file for further reference.
//

#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserParallelWorld.hh"
#include "sstream"
#include "G4Transform3D.hh"


#ifndef PMTReadoutGeometry_h
#define PMTReadoutGeometry_h


namespace larg4 
{

  class PMTReadoutGeometry : public G4VUserParallelWorld
    {
    public:
      PMTReadoutGeometry(G4String PMTSensitiveName, const G4String name = "PMTReadoutGeometry");
      virtual ~PMTReadoutGeometry();
      
      virtual void Construct();
    private:
      void                            FindVolumes(G4VPhysicalVolume *, G4String, std::vector<G4Transform3D>);
      std::vector<G4LogicalVolume*>   PMTVolumes;
      std::vector<G4Transform3D>      PMTTransformations;
      G4String fPMTSensitiveName;

    };
  
}

#endif
