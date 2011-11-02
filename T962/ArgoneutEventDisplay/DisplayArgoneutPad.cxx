///
/// \file    DisplayArgoneutPad.cxx
/// \brief   Drawing pad showing a 3D rendering of the ArgoNeuT detector
/// \author  msoderbe@syr.edu
///
#include "T962/ArgoneutEventDisplay/DisplayArgoneutPad.h"
#include <iostream>
#include "TPad.h"
#include "TView3D.h"
#include "TGLViewer.h"
#include "EventDisplayBase/evdb.h"
#include "EventDisplayBase/EventHolder.h"
#include "Geometry/Geometry.h"
#include "EventDisplay/GeometryDrawer.h"
#include "EventDisplay/SimulationDrawer.h"
#include "T962/ArgoneutEventDisplay/ArgoneutRecoBaseDrawer.h"
#include "T962/ArgoneutEventDisplay/MinosDrawer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

using namespace evd;
using namespace argoevd;

///
/// Create a pad to show a 3D rendering of the detector and events
/// @param nm : Name of the pad
/// @param ti : Title of the pad
/// @param x1 : Location of left  edge of pad (0-1)
/// @param x2 : Location of right edge of pad (0-1)
/// @param y1 : Location of bottom edge of pad (0-1)
/// @param y2 : Location of top    edge of pad (0-1)
/// @param opt: Options. Currently just a place holder
///
DisplayArgoneutPad::DisplayArgoneutPad(const char* nm, const char* ti,
                                       double x1, double y1,
                                       double x2, double y2) : 
   DrawingPad(nm, ti, x1, y1, x2, y2)
{
  this->Pad()->SetFillColor(kBlack);
  this->Pad()->Draw();
  this->Pad()->cd();
  fView = new evdb::View3D();
}

//......................................................................

DisplayArgoneutPad::~DisplayArgoneutPad() 
{
  if (fView) { delete fView; fView = 0; }
  if (fArgoneutRecoDraw){delete fArgoneutRecoDraw; fArgoneutRecoDraw=0;}
  if (fMinosDraw){delete fMinosDraw; fMinosDraw=0;}
     
}

//......................................................................

void DisplayArgoneutPad::Draw() 
{
  fView->Clear();

  art::ServiceHandle<geo::Geometry> geo;

  // grab the event from the singleton
  const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();

  if(evt){
    this->ArgoneutRecoBaseDraw()->  ArgoProng3D     (*evt, fView);
    this->MinosDraw()->  DetOutline3D(fView);
    this->MinosDraw()->  Minos3D     (*evt, fView);
    this->MinosDraw()->  MinosOutline3D(fView);
  }
 

  this->Pad()->Clear();
  this->Pad()->cd();
  if (fPad->GetView()==0) {
    int irep;
    double rmin[]={-2.1*geo->DetHalfWidth(),-2.1*geo->DetHalfHeight(),-0.5*geo->DetLength()};
    double rmax[]={ 2.1*geo->DetHalfWidth(), 2.1*geo->DetHalfHeight(), 2.1*geo->DetLength()};
    TView3D* v = new TView3D(1,rmin,rmax);
    v->SetPerspective();
    v->SetView(0.0,260.0,270.0,irep);
    fPad->SetView(v); // ROOT takes ownership of object *v
  }
  fView->Draw();
  fPad->Update();
}

ArgoneutRecoBaseDrawer* DisplayArgoneutPad::ArgoneutRecoBaseDraw() 
{
   fArgoneutRecoDraw = new ArgoneutRecoBaseDrawer();
   return fArgoneutRecoDraw;
}

MinosDrawer* DisplayArgoneutPad::MinosDraw() 
{
   fMinosDraw = new MinosDrawer();
   return fMinosDraw;
}



////////////////////////////////////////////////////////////////////////
