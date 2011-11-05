////////////////////////////////////////////////////////////////////////
/// \file PMTLookup.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Provide a map between G4VPhysicalVolumes of PMTs and PMT ID's.
//
// At the beginning of the event, DetectorConstruction passes us
// a pointer to the physical volume store.  The SetPhysicalVolumes
// method loops through, finding the PMT volumes.
// 
// Initially these are all named identically by the GEANT4 gdml parser.
// We rename them so that each is unique, and use the physical vol names 
// as references to the individual PMTs. A map between PMT ID and volume
// name is created.
//
// We can then query this map using the GetID function, querying either
// by physical volume or by name.
//
// The GetN() function gives the number of PMT volumes which have been
// registered.
// 
// Another option would be to use the addresses of the volumes in the pvs
// as the key.  This would be perhaps slightly faster, but  
// somewhat less transparent.  I may change to this scheme later if
// it will significantly help performance.
//
// Ben Jones, MIT, 06/04/2010
//

#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include <map>

#ifndef PMTLOOKUP_h
#define PMTLOOKUP_h 1


namespace larg4 {
  class PMTLookup
    {
    public:
      ~PMTLookup(){}
      static PMTLookup * Instance();
      void SetPhysicalVolumes(G4PhysicalVolumeStore *, std::string PMTVolLabel);
      void AddPhysicalVolume(G4VPhysicalVolume *);
      int GetID(G4VPhysicalVolume *);
      int GetID(std::string);
      int GetN();
      
    protected:
      PMTLookup();

    private:
      std::map<std::string, int> fTheIDMap;
      int fTheTopID;
      
    };

}


#endif
