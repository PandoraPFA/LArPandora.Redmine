////////////////////////////////////////////////////////////////////////
/// \file    PhotonLibraryObjects.h
/// \version $Id: CalHit.cxx,v 1.9 2011/01/24 23:18:19 p-novaart Exp $
/// \author  bpjones
////////////////////////////////////////////////////////////////////////
#ifndef PhotonLibraryObjects_h
#define PhotonLibraryObjects_h 1

#include "Simulation/PMTHit.h"
#include "Simulation/PhotonVoxels.h"
#include "TObject.h"
#include "Utilities/VectorMap.h"

namespace phot {
  
  class TPMTHitVMap : public TObject, public util::VectorMap<int, sim::PMTHitCollection>
  {
  public:
    TPMTHitVMap();
    virtual ~TPMTHitVMap();
    
    ClassDef(TPMTHitVMap,1);
  };
  
  class TPhotonVoxelDef : public TObject, public sim::PhotonVoxelDef
  {
  public:
    TPhotonVoxelDef();
    TPhotonVoxelDef(sim::PhotonVoxelDef);
    virtual ~TPhotonVoxelDef();
    
    ClassDef(TPhotonVoxelDef,1);
  };

  class TSampleSizeV : public TObject, public std::vector<Int_t>
  {
  public:
    TSampleSizeV();
    virtual ~TSampleSizeV();
    
    ClassDef(TSampleSizeV,1);

  };


}


#endif

#ifdef __CINT__
#include "Simulation/PhotonVoxels.h"
#include "Simulation/PMTHit.h"
#endif
