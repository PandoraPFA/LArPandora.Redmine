////////////////////////////////////////////////////////////////////////
/// \file    PhotonLibraryObjects.cxx
/// \version $Id: CalHit.cxx,v 1.9 2011/01/24 23:18:19 p-novaart Exp $
/// \author  bpjones
////////////////////////////////////////////////////////////////////////

#include "PhotonPropagation/PhotonLibraryObjects.h"

ClassImp(phot::TPMTHitVMap);
ClassImp(phot::TPhotonVoxelDef);
ClassImp(phot::TSampleSizeV);

namespace phot {


  TPMTHitVMap::TPMTHitVMap()
  {
  }
  
  TPMTHitVMap::~TPMTHitVMap()
  {
  }
  
  TPhotonVoxelDef::TPhotonVoxelDef():sim::PhotonVoxelDef()
  {
  }
  
  TPhotonVoxelDef::~TPhotonVoxelDef()
  {
  }
  
  TPhotonVoxelDef::TPhotonVoxelDef(sim::PhotonVoxelDef pvd): sim::PhotonVoxelDef(pvd)
  {
  }
  
  TSampleSizeV::TSampleSizeV()
  {
  }

  TSampleSizeV::~TSampleSizeV()
  {
  }



}
