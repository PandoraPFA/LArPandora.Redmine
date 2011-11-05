////////////////////////////////////////////////////////////////////////
/// \file    PhotonLibrary.h
/// \version $Id: CalHit.cxx,v 1.9 2011/01/24 23:18:19 p-novaart Exp $
/// \author  bpjones
////////////////////////////////////////////////////////////////////////
#include "Simulation/PhotonVoxels.h"
#include "Simulation/PMTHit.h"
#include "Utilities/VectorMap.h"
#include "TVector3.h"
#include "PhotonLibraryObjects.h"
#include "TFile.h"

#ifndef PhotonLibrary_h
#define PhotonLibrary_h 1

namespace phot {

  class PhotonLibrary : public TObject
  {
  public:
   
    PhotonLibrary();
    PhotonLibrary(sim::PhotonVoxelDef VoxelDef);

    virtual ~PhotonLibrary();
    
    // Geometry Accessor
    virtual sim::PhotonVoxelDef GetVoxelDef();

    // PMTHit Accessors
    sim::PMTHitCollection  GetHitCollectionForVoxel( int ID );
    sim::PMTHitCollection  GetHitCollectionForVoxel( TVector3 Position );

    // Setters and Incrementers
    void SetHitCollectionForVoxel  (int ID, sim::PMTHitCollection, int);
    void AddToHitCollectionForVoxel(int ID, sim::PMTHitCollection, int );

    void SaveInFile(TFile * TF);
    void LoadFromFile(TFile * TF);

  private:
    TPMTHitVMap                  * fHitCollectionsForVoxels; //->
    TPhotonVoxelDef              * fPhotonVoxelDef; //->
    TSampleSizeV                 * fSampleSizesForVoxels; //->
    
    ClassDef(PhotonLibrary,4); 
  };
}

#endif
