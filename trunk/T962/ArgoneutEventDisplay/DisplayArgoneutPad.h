///
/// \file    DisplayArgoneutPad.h
/// \brief   Drawing pad showing a 3D rendering of the ArgoNeuT detector
/// \author  msoderbe@syr.edu
///
#ifndef EVD_DISPLAYARGONEUTPAD_H
#define EVD_DISPLAYARGONEUTPAD_H
#include "EventDisplay/DrawingPad.h"
class TH3F;
namespace evdb { class View3D; }

namespace argoevd {
   class ArgoneutRecoBaseDrawer;
   class MinosDrawer;
  
  /// A drawing pad showing a 3D rendering of the detector
   class DisplayArgoneutPad : public evd::DrawingPad {
  public:
      DisplayArgoneutPad(const char* nm, const char* ti,
                         double x1, double y1,
                         double x2, double y2);
      ~DisplayArgoneutPad();
      void Draw();
      
      ArgoneutRecoBaseDrawer* ArgoneutRecoBaseDraw();
      MinosDrawer*  MinosDraw();  
      
  private:
    // void AutoZoom();
  private:
    evdb::View3D* fView;  ///< Collection of graphics objects to render
    ArgoneutRecoBaseDrawer* fArgoneutRecoDraw;   ///< Drawer for recobase objects  
    MinosDrawer*            fMinosDraw;   ///< Drawer for recobase objects      
      
  };
};

#endif
////////////////////////////////////////////////////////////////////////
