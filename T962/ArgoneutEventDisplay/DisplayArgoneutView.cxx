///
/// \file    DisplayArgoneutView.cxx
/// \brief   The "main" event display view that most people will want to use
/// \author  msoderbe@syr.edu
///
#include "TCanvas.h"
#include "TVirtualViewer3D.h"
#include "T962/ArgoneutEventDisplay/DisplayArgoneutView.h"
#include "T962/ArgoneutEventDisplay/DisplayArgoneutPad.h"

#include "art/Framework/Principal/Event.h"

namespace argoevd{

  //......................................................................
  DisplayArgoneutView::DisplayArgoneutView(TGMainFrame* mf) : evdb::Canvas(mf)
  {
    evdb::Canvas::fCanvas->cd();
    
    fDisplayArgoneutPad = new DisplayArgoneutPad("fDisplayArgoneutPad","3D Display",
				     0.0, 0.0, 1.0, 1.0);
    
    this->Connect("CloseWindow()","argoevd::DisplayArgoneutView",this,"CloseWindow()");
    
    //fDisplayArgoneutPad->Draw();
    
    evdb::Canvas::fCanvas->Update();
  }
  
  //......................................................................
  DisplayArgoneutView::~DisplayArgoneutView() 
  {
  }

  //......................................................................
  void DisplayArgoneutView::CloseWindow()
  {
    delete this;
  }

  //......................................................................
  void DisplayArgoneutView::Draw(const char* /*opt*/) 
  {
    fDisplayArgoneutPad->Draw();
    evdb::Canvas::fCanvas->Update();
    
    TVirtualViewer3D *viewer = fDisplayArgoneutPad->Pad()->GetViewer3D("ogl");
    //TVirtualViewer3D *viewer = fDisplayArgoneutPad->Pad()->GetViewer3D("x3d");
    viewer->PreferLocalFrame();
    viewer->ResetCameras();
    viewer->PadPaint(fDisplayArgoneutPad->Pad());
    
  }

}// end namespace
////////////////////////////////////////////////////////////////////////
