#include "iostream"
#include "PhotonLibraryParameters.h"


namespace sim {
  PhotonLibraryParameters::PhotonLibraryParameters()
  {
  }

  PhotonLibraryParameters::PhotonLibraryParameters(PhotonVoxelDef  theDef, int theID, int theCount)
  {
    fThePhotonVoxelDefinition  =  theDef;
    fTheVoxelID                =  theID;
    fThePhotonCount            =  theCount;
  }

  PhotonLibraryParameters::~PhotonLibraryParameters()
  {
  }

}
