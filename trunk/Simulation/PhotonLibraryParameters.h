#include "PhotonVoxels.h"
#include "TObject.h"

#ifndef PhotonLibraryParameters_h
#define PhotonLibraryParameters_h 1

namespace sim {
  
  class PhotonLibraryParameters
    {
    public:

      PhotonLibraryParameters();
      PhotonLibraryParameters(PhotonVoxelDef, int, int);
      virtual ~PhotonLibraryParameters();
      
      PhotonVoxelDef   GetVoxelDef() const {return fThePhotonVoxelDefinition;}
      int            GetVoxelID()        const {return fTheVoxelID;}
      int            GetPhotonCount()    const {return fThePhotonCount;}
      
    private:
      PhotonVoxelDef      fThePhotonVoxelDefinition;   
      int               fTheVoxelID;                 
      int               fThePhotonCount;            


    };


}

#endif
