//////////////////////////////////////////////////////////////////////////
/// \file MagField.cxx
///
/// \version $Id$
/// \author dmckee@phys.ksu.edu
//////////////////////////////////////////////////////////////////////////
/// \class MagField MagField.h 
/// The initial implementation will be trivial: simply supporting a
/// constant field in a named detector volume. In principle we should
/// read a full field map from an external file of some kind.
///
/// We support three FHICL value for now:
///
///    - "UseField" a boolean. When false we don't even instantiate a
///      Magnetic field object
///    - "Constant Field" a vector< double > which should have three
///      elements and is interpreted in Tesla
///    - "MagnetizedVolume" names the G4logical volume to which the
///      field should be attached
//////////////////////////////////////////////////////////////////////////

#include "MagField/mag.h"

#include <vector>
#include <string>

namespace mag {

  MagField::MagField(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
    :fUseField(pset.get< bool >("UseField"))
  {
    // These need to be read as types that FHICL know about, but they
    // are used by Geant, so I store them in Geant4 types.
    std::vector<double> field  = 
      (pset.get<std::vector<double> >("ConstantField"));
    std::string         volume = 
      (pset.get<std::string         >("MagnetizedVolume"));
    // Force the dimension of the feild definition
    field.resize(3);
    for (int i=0; i<3; i++) fField[i] = field[i];
    fVolume = volume.c_str();
  }

  G4ThreeVector MagField::FieldAtPoint(G4ThreeVector p)const{
    // FIXME: This does not do what it says. Must test to see if the
    // point is in the master volume
    //
    // But it is enough to let me code the DetectorConstruction bit

    if ( /* is in the magnetized volume */ true ) return fField;
    return G4ThreeVector(0);
  }
  G4String MagField::MagnetizedVolume()const{
    return fVolume;
  }

}
