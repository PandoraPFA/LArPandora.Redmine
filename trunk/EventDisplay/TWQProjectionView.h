////////////////////////////////////////////////////////////////////////
///
/// \file    TWQProjectionView.h
/// \brief   A view showing the time vs wire, charge and charge vs time information for an event
/// \author  brebel@fnal.gov
/// \version $Id: TQPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
////////////////////////////////////////////////////////////////////////

#ifndef EVD_TWQPROJECTION_H
#define EVD_TWQPROJECTION_H

#include "EventDisplayBase/Canvas.h"
#include "RQ_OBJECT.h"
#include <vector>

// Forward declarations
class TGCheckButton;
class TGRadioButton;
class TGNumberEntry;
class TGLabel;
class TRootEmbeddedCanvas;

namespace evd {

  class MCBriefPad;
  class TQPad;
  class TWireProjPad;
  class HeaderPad;

  class TWQProjectionView : public evdb::Canvas {

  public:
    
    RQ_OBJECT("evd::TWQProjectionView");

  public:
    TWQProjectionView(TGMainFrame* mf);
    ~TWQProjectionView();

    const char* Description() const { return "Time/Wire/Charge Projections"; }
    const char* PrintTag()    const { return "twq-proj";                     }

    void Draw(const char* opt="");

 //   void        RangeChanged();
    static void ChangeWire(int plane, void * wqpv); 
    static void SetMouseZoomRegion(int plane, void * wqpv);
    void        SetPlaneWire();
    void        SetPlane();
    void        SetWire();
    void        SetThreshold();
    void        SetGreyscale();
    void        SetMCInfo();
    void        SetRawCalib();
    void 	SetUpSideBar();
    void        SetUpZoomButtons();
    void        SetZoom(int plane,int wirelow,int wirehi,int timelo,int timehi);
    void 	ZoomInterest(bool flag=true);
    void 	SetZoomInterest();
    void 	PrintCharge();
    
    
    
  private:

    HeaderPad*  fHeaderPad;              ///< Show header information    
    TQPad*      fWireQ;                  ///< Histogram of charge vs time on selected wire
    MCBriefPad* fMC;                     ///< Short summary of MC event    
    std::vector<TQPad* >        fPlaneQ; ///< charge on each plane
    std::vector<TWireProjPad*>  fPlanes; ///< time vs wire projection for each plane

    TGCompositeFrame*    fVFrame;        ///< needed for the side frame
    TGCompositeFrame*    fMetaFrame;   ///< needed for the side frame

    TGLabel* fWireLabel;
    TGLabel* fPlaneLabel;
    TGLabel* fThresLabel;
    TGLabel* fGreyLabel;

    TGNumberEntry* fWireEntry;     ///< Wire number displayed.	
    TGNumberEntry* fPlaneEntry;    ///< Plane number displayed.	
    TGNumberEntry* fThresEntry;    ///< ADC threshold to display.	
    TGCheckButton* fGreyScale;     ///< Display gray or color scale
    TGCheckButton* fMCOn;          ///< Display MC truth information
    TGRadioButton* fRawDraw;       ///< Draw Raw information only
    TGRadioButton* fCalibDraw;     ///< Draw calibrated information only
    TGRadioButton* fRawCalibDraw;  ///< Draw raw and calibrated information

   

    TGTextButton* fZoomInterest;       ///< Zoom on iteresting region
    TGTextButton* fUnZoomInterest;       ///< Unzoom on iteresting region
    TGCheckButton* fToggleAutoZoom;       ///< Toggle the autozoom setting 
  };

}// namespace

#endif //EVD_TWQPROJECTION_H
