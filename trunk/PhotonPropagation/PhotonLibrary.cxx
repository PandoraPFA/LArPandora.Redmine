////////////////////////////////////////////////////////////////////////
/// \file    PhotonLibrary.cxx
/// \version $Id: CalHit.cxx,v 1.9 2011/01/24 23:18:19 p-novaart Exp $
/// \author  bpjones
////////////////////////////////////////////////////////////////////////
#include "PhotonPropagation/PhotonLibrary.h"
#include "Simulation/PhotonVoxels.h"

namespace phot{

  PhotonLibrary::PhotonLibrary()
  {
    std::cout<<"PhotonLibrary Default constructor, to keep root happy" <<std::endl;
    fHitCollectionsForVoxels = new TPMTHitVMap();
    fPhotonVoxelDef = new TPhotonVoxelDef();
  }

  PhotonLibrary::PhotonLibrary(sim::PhotonVoxelDef VoxelDef)
  {
    fPhotonVoxelDef = new TPhotonVoxelDef( VoxelDef);
    fHitCollectionsForVoxels = new TPMTHitVMap();
    fSampleSizesForVoxels = new TSampleSizeV();
  }

  PhotonLibrary::~PhotonLibrary()
  {
  }

  sim::PhotonVoxelDef PhotonLibrary::GetVoxelDef()
  {
    return *fPhotonVoxelDef;
  }

  void PhotonLibrary::SaveInFile(TFile * TF)
  {
    TF->WriteObject(fHitCollectionsForVoxels,"TheHitCollection");
    TF->WriteObject(fPhotonVoxelDef,"TheVoxelDef");      
    TF->WriteObject(fSampleSizesForVoxels,"TheSampleSizes");
  }

  void PhotonLibrary::LoadFromFile(TFile * TF)
  {
    fPhotonVoxelDef          = (TPhotonVoxelDef*) TF->Get("TheVoxelDef");
    fHitCollectionsForVoxels = (TPMTHitVMap*)     TF->Get("TheHitCollection");
    fSampleSizesForVoxels     = (TSampleSizeV*)    TF->Get("TheSampleSizes"); 
  }

  sim::PMTHitCollection PhotonLibrary::GetHitCollectionForVoxel(int ID)
  {
    if(fPhotonVoxelDef->IsLegalVoxelID(ID))
      return (*fHitCollectionsForVoxels)[ID];
    else
      {
	std::cout<<"PhotonLibrary Error - voxel ID out of range" <<std::endl;
	throw("");
      }
  }

  sim::PMTHitCollection PhotonLibrary::GetHitCollectionForVoxel(TVector3 Position)
  {
    return GetHitCollectionForVoxel(fPhotonVoxelDef->GetVoxelID(Position));
  }
  
  void PhotonLibrary::SetHitCollectionForVoxel(int ID,sim::PMTHitCollection TheCollection, int SampleSize)
  {
    (*fHitCollectionsForVoxels)[ID] = TheCollection;
    (*fSampleSizesForVoxels)[ID] = SampleSize; 
  }

  void PhotonLibrary::AddToHitCollectionForVoxel(int ID,sim::PMTHitCollection TheCollection, int SampleSize)
  {
    ((*fHitCollectionsForVoxels)[ID]) += TheCollection;
    (*fSampleSizesForVoxels)[ID] += SampleSize;
  }

}
