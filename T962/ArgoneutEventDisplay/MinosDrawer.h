#ifndef ARGONEUTEVD_MINOSDRAWER_H
#define ARGONEUTEVD_MINOSDRAWER_H
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

namespace evdb{ 
  class View2D;
  class View3D;
}
namespace geo      { class Geometry; }
namespace t962     { class MINOS; } 



namespace argoevd {
  /// Aid in the rendering of Geometry objects
  class MinosDrawer {
  public:
    MinosDrawer();
    ~MinosDrawer();
    void DetOutline3D(evdb::View3D* view);
     void MinosOutline3D(evdb::View3D* view);

     void Minos3D(const art::Event& evt,
                  evdb::View3D*     view);
     void MinosTrack(const art::Ptr<t962::MINOS> minos,
                     int                 id,
                     evdb::View3D*       view,
                     bool                matched);

  private:
     int GetMinos(const art::Event&            evt,
                  const std::string&           which,
                  art::PtrVector<t962::MINOS>& minos);

     
  };
};

#endif
////////////////////////////////////////////////////////////////////////
