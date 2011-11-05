#ifndef __CINT__

#include <TString.h>
#include <TVector3.h>
#include <Rtypes.h>
#include <TFile.h>

#include <vector>
#include <map>
#include <cstring>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Simulation/PMTHit.h"
#include "Simulation/PhotonVoxels.h"
#include "PhotonPropagation/PhotonLibrary.h"

#ifndef PhotonLibraryService_h
#define PhotonLibraryService_h

namespace phot
{

  class PhotonLibraryService
  {
  public:
    PhotonLibraryService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~PhotonLibraryService();

    PhotonLibrary* GetLibrary();
    
    void CreateLibraryFile(std::string, sim::PhotonVoxelDef);
    void OpenLibraryFile(std::string);
    
    void CheckVoxelDef(sim::PhotonVoxelDef);

    void InitializeBuildJob(sim::PhotonVoxelDef, bool ExtendExisting=false);
    void InitializeReadJob();

  private:
    // Parameters from config
    std::string fFileName;
    bool fExtendExisting;
    bool fBuildLibrary;
  
    //Voxel geometry definition
    sim::PhotonVoxelDef fVoxelDef;
   
    //Pointer to library file
    TFile *fTheLibraryFile;
    
    //Pointer to the library
    PhotonLibrary * fTheLibrary;
  };



}

#endif

#endif
