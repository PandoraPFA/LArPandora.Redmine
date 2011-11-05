/// \file MagField.h
/// \brief Describe the magnetic field structure of a detector
/// 
/// \version $Id$
/// \author dmckee@phys.ksu.edu
//////////////////////////////////////////////////////////////////////////
/// \namespace mag
/// A namespace for simulated magnetic fields
//////////////////////////////////////////////////////////////////////////
#ifndef MAG_FIELD_H
#define MAG_FIELD_H

// Geant4 includes
#include <G4String.hh>
#include <G4ThreeVector.hh>

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace mag {

  // Specifies the magnetic field over all space
  //
  // The default implementation, however, uses a nearly trivial,
  // non-physical hack.
  class MagField {
  public:
    MagField(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~MagField(){};

    bool UseField()const{return fUseField;};
    // return the field at a particular point
    G4ThreeVector FieldAtPoint(G4ThreeVector p=G4ThreeVector(0))const;
    // return the outermost affected volume
    G4String MagnetizedVolume()const;

  private:
    // The simplest implmentation has a constant field inside a named
    // detector volume
    bool fUseField;
    G4ThreeVector fField;
    G4String fVolume;
  };

}

#endif // MAG_FIELD_H
