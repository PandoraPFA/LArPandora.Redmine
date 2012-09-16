#ifndef ARGONEUTEVD_PADDLESDRAWER_H
#define ARGONEUTEVD_PADDLESDRAWER_H
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
namespace t962     { class Paddles; } 



namespace argoevd {
  
  class PaddlesDrawer {
  public:
    PaddlesDrawer();
    ~PaddlesDrawer();
    
    void PaddlesOutline3D(evdb::View3D* view,
                          int paddlenumber,
                          int color);

    void DrawPaddlesInfo3D(const art::Event& evt,
                       evdb::View3D*     view);

    void FindCoincidences(const art::Handle< t962::Paddles > p, int* color, int tol);
 
  private:
 

     
  };
};

#endif
////////////////////////////////////////////////////////////////////////
