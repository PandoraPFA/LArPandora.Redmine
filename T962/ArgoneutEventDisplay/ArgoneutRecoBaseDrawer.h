/// \file    ArgoNeuTRecoBaseDrawer.h
/// \brief   Class to aid in the rendering of ArgoNeuT Reconstruction objects
/// \author  msoderbe@syr.edu
#ifndef EVD_ARGONEUTRECOBASEDRAWER_H
#define EVD_ARGONEUTRECOBASEDRAWER_H
#include "EventDisplay/RecoBaseDrawer.h"
#include <vector>
#ifdef __CINT__
namespace art { 
  class Event;
  class PtrVector;
  class Ptr;
}
#else
#include "art/Persistency/Common/PtrVector.h"
#endif

namespace t962     { class MINOS; }

namespace argoevd {
  /// Aid in the rendering of ArgoNeuT Reconstruction objects
   class ArgoneutRecoBaseDrawer : public evd::RecoBaseDrawer {
  public:
    ArgoneutRecoBaseDrawer();
    ~ArgoneutRecoBaseDrawer();

  public:

    void ArgoSpacePoint(std::vector<const recob::SpacePoint*> sps,
                         int                 id,
                         evdb::View3D*       view,
                         bool matched=false,
                         float charge=-999.);
     
     void ArgoProng3D(const art::Event& evt,
                      evdb::View3D*     view);

    
  private:

    int GetArgoTracks(const art::Event&                 evt,
                      const std::string&                which,
                      art::PtrVector<recob::Track>&     track);
    int GetArgoShowers(const art::Event&                evt,
		       const std::string&               which,
		       art::PtrVector<recob::Shower>&   shower);


     int GetMinos(const art::Event&            evt,
                  const std::string&           which,
                  art::PtrVector<t962::MINOS>& minos);


    
        
    
    
  };
}

#endif
////////////////////////////////////////////////////////////////////////
