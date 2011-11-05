///
/// \file    Display3DView.cxx
/// \brief   The "main" event display view that most people will want to use
/// \author  messier@indiana.edu
/// \version $Id: Display3DView.cxx,v 1.2 2011/01/19 16:40:59 p-novaart Exp $
///
#include "TCanvas.h"
#include "TVirtualViewer3D.h"
#include "EventDisplay/Display3DView.h"
#include "EventDisplay/Display3DPad.h"

#include "art/Framework/Principal/Event.h"

namespace evd{

  //......................................................................
  Display3DView::Display3DView(TGMainFrame* mf) : evdb::Canvas(mf)
  {
    evdb::Canvas::fCanvas->cd();
    
    fDisplay3DPad = new Display3DPad("fDisplay3DPad","3D Display",
				     0.0, 0.0, 1.0, 1.0, "");
    
    this->Connect("CloseWindow()","evd::Display3DView",this,"CloseWindow()");
    
    fDisplay3DPad->Draw();
    
    evdb::Canvas::fCanvas->Update();
  }
  
  //......................................................................
  Display3DView::~Display3DView() 
  {
  }

  //......................................................................
  void Display3DView::CloseWindow()
  {
    delete this;
  }

  //......................................................................
  void Display3DView::Draw(const char* /*opt*/) 
  {
    fDisplay3DPad->Draw();
    evdb::Canvas::fCanvas->Update();
    
    TVirtualViewer3D *viewer = fDisplay3DPad->Pad()->GetViewer3D("ogl");
    viewer->PreferLocalFrame();
    viewer->ResetCameras();
    viewer->PadPaint(fDisplay3DPad->Pad());
    
  }

}// end namespace
////////////////////////////////////////////////////////////////////////
