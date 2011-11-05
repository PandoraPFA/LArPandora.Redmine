/////////////////////////////////////////////////////////////////////////////
///
/// \file    TWireProjPad.h
/// \brief   Drawing pad showing a single X-Z or Y-Z projection of an event
/// \author  messier@indiana.edu
/// \version $Id: TWireProjPad.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
///
/////////////////////////////////////////////////////////////////////////////

#ifndef EVD_TWIREPROJPAD_H
#define EVD_TWIREPROJPAD_H
#include "EventDisplay/DrawingPad.h"

class TH1F;

namespace evdb { class View2D;   }

namespace evd {

  /// A drawing pad for time vs wire
  class TWireProjPad : public DrawingPad {
  public:
    TWireProjPad(const char* nm, const char* ti,
		 double x1, double y1,
		 double x2, double y2,
		 unsigned int plane);
    ~TWireProjPad();
    void Draw(const char* opt=0);
    
    void GetWireRange(int *i1, int *i2) const;
    void SetWireRange(int i1, int i2);

    void SetZoomRange(int i1, int i2,int y1, int y2);
    
    unsigned int GetPlane() const { return fPlane; }

     void ShowFull(int override=0);
  private:
/*     void AutoZoom(); */
   

  private:

    unsigned int  fPlane; ///< Which plane in the detector
    TH1F*         fHisto; ///< Histogram to draw object on
    evdb::View2D* fView;  ///< Collection of graphics objects to render

    double        fXLo;   ///< Low  value of x axis
    double        fXHi;   ///< High value of x axis
    double        fYLo;   ///< Low  value of y axis
    double        fYHi;   ///< High value of y axis
    int           fOri;   ///< Orientation of the axes - see RawDrawingOptions for values
  };
};

#endif
////////////////////////////////////////////////////////////////////////
