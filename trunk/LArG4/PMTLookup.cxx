////////////////////////////////////////////////////////////////////////
/// \file PMTLookup.cxx
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the PMTLookup class.
//
// See comments in the PMTLookup.h file.
//
// Ben Jones, MIT, 06/04/2010
//


#include "LArG4/PMTLookup.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {
  PMTLookup * ThePMTLookup;
  
  //--------------------------------------------------
  PMTLookup::PMTLookup()
  {
    fTheTopID=0;
  }

  //--------------------------------------------------
  PMTLookup * PMTLookup::Instance()
  {
    if(!ThePMTLookup){
      ThePMTLookup = new PMTLookup;
    }
    return ThePMTLookup;  
  }

  //--------------------------------------------------
  int PMTLookup::GetID(std::string TheName)
  {
    int ID = fTheIDMap[TheName];
    if(ID==0){
      throw cet::exception("PMTLookup") <<"PMTLookup : Attempted to lookup a PMT "
					<< "with unknown name : " 
					<< TheName.c_str();
    }
    else
      return ID;
  }

  //--------------------------------------------------
  int PMTLookup::GetID(G4VPhysicalVolume* TheVolume)
  {
    std::string TheName = TheVolume->GetName();
    return GetID(TheName);
  }

  //--------------------------------------------------
  void PMTLookup::AddPhysicalVolume(G4VPhysicalVolume * volume)
  {
    // If volume with this name already exists, rename it
    fTheTopID++;
    G4String VolName = volume->GetName();
    if(fTheIDMap[VolName]){
      std::stringstream VolName;
      VolName << volume->GetName() << "_" << fTheTopID;
      volume->SetName(VolName.str().c_str());
      fTheIDMap[VolName.str()] = fTheTopID;
    }
    // Else just add it to the lookup table
    else{
      fTheIDMap[VolName] = fTheTopID;
    }
  }

  //--------------------------------------------------
  void PMTLookup::SetPhysicalVolumes(G4PhysicalVolumeStore * pvs, std::string PMTVolLabel)
  {
    for ( G4PhysicalVolumeStore::iterator i = pvs->begin(); i != pvs->end(); ++i ){
      G4VPhysicalVolume* volume = (*i);
      G4String name = volume->GetName();
      if(name.contains(PMTVolLabel.c_str())){
	fTheTopID++;
	std::stringstream VolName;
	VolName << volume->GetName() << "_" << fTheTopID;
	volume->SetName(VolName.str().c_str());
	fTheIDMap[VolName.str()] = fTheTopID;
      }
    }

  }

  //--------------------------------------------------
  int PMTLookup::GetN()
  {
    return fTheTopID;
  }
  
}
