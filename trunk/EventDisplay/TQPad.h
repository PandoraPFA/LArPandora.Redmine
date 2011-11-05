///
/// \file    TQPad.h
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
/// \version $Id: TQPad.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
///
#ifndef EVD_TQPAD_H
#define EVD_TQPAD_H
#include "EventDisplay/DrawingPad.h"
namespace evdb { class View2D; }
class TH1F;

namespace evd {
  class TQPad : public DrawingPad{
  public:
    TQPad(const char* nm, const char* ti,
	  double x1, double y1,
	  double x2, double y2,
	  const char *opt,
	  unsigned int plane,
	  unsigned int wire);
    ~TQPad();
    void Draw();
    
    void SetPlaneWire(unsigned int plane=0, unsigned int wire=0) { fPlane = plane; fWire = wire; }

  private:
    void BookHistogram();

    unsigned int  fWire;
    unsigned int  fPlane;     ///< Which plane in the detector
    int           fTQ;        ///< 0 = plot shows charge only, 1 = plot shows charge vs time for a wire
    TH1F*         fRawHisto;  ///< 1-D Histogram of charge or charge vs time
    TH1F*         fRecoHisto; ///< 1-D Histogram of charge or charge vs time
    evdb::View2D* fView;      ///< Superimpose scale on 1D histo
  };
}

#endif
////////////////////////////////////////////////////////////////////////
