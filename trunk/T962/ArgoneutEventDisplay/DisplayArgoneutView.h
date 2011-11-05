///
/// \file    DisplayArgoneutView.h
/// \brief   A view showing a 3D rendering of the ArgoNeuT detector
/// \author  msoderbe@syr.edu
///
#ifndef EVD_DISPLAYARGONEUTVIEW_H
#define EVD_DISPLAYARGONEUTVIEW_H
#include "EventDisplayBase/evdb.h"
#include "RQ_OBJECT.h"

namespace argoevd {
  class DisplayArgoneutPad;

  /// View of event shoing the XZ and YZ readout planes
  class DisplayArgoneutView : public evdb::Canvas {

  public:
    
    RQ_OBJECT("argoevd::DisplayArgoneutView");

  public:
    DisplayArgoneutView(TGMainFrame* mf);
    ~DisplayArgoneutView();
    
    const char* Description() const { return "3D ArgoNeuT Display"; }
    const char* PrintTag()    const { return "larargoneut";               }
    void Draw(const char* opt="");
    void CloseWindow();

  private:
    DisplayArgoneutPad* fDisplayArgoneutPad; /// Pad showing 3D view of the detector
  };
}

#endif
////////////////////////////////////////////////////////////////////////
