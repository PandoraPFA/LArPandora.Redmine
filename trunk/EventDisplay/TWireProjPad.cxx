////////////////////////////////////////////////////////////////////////
///
/// \file    TWireProjPad.cxx
/// \brief   Drawing pad for X-Z or Y-Z projections of events
/// \author  messier@indiana.edu
/// \version $Id: TZProjPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
////////////////////////////////////////////////////////////////////////
#include <algorithm>

#include "EventDisplay/TWireProjPad.h"
#include "TPad.h"
#include "TH1F.h"
#include "TString.h"
#include "EventDisplayBase/evdb.h"
#include "EventDisplayBase/EventHolder.h"
#include "Geometry/geo.h"
#include "EventDisplay/GeometryDrawer.h"
#include "EventDisplay/SimulationDrawer.h"
#include "EventDisplay/RawDataDrawer.h"
#include "EventDisplay/RecoBaseDrawer.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "EventDisplay/EvdLayoutOptions.h"
#include "TFrame.h"


#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"

namespace evd{

  ///
  /// Create a pad showing a single X-Z or Y-Z projection of the detector
  /// \param nm : Name of the pad
  /// \param ti : Title of the pad
  /// \param x1 : Location of left  edge of pad (0-1)
  /// \param x2 : Location of right edge of pad (0-1)
  /// \param y1 : Location of bottom edge of pad (0-1)
  /// \param y2 : Location of top    edge of pad (0-1)
  /// \param plane : plane number of view
  ///
  TWireProjPad::TWireProjPad(const char* nm,
			     const char* ti,
			     double x1, double x2,
			     double y1, double y2,
			     unsigned int plane)
    : DrawingPad(nm, ti, x1, x2, y1, y2)
    , fPlane(plane)
  {
  
    art::ServiceHandle<geo::Geometry> geo;

    this->Pad()->SetBit(kCannotPick);
    this->Pad()->SetBit(TPad::kCannotMove);
    this->Pad()->cd();

    this->Pad()->SetLeftMargin  (0.070);
    this->Pad()->SetRightMargin (0.010);

    // how many planes in the detector and 
    // which plane is this one?
    
    unsigned int planes = geo->Nplanes();
    this->Pad()->SetTopMargin   (0.005);
    this->Pad()->SetBottomMargin(0.110);

    // there has to be a better way of doing this that does
    // not have a case for each number of planes in a detector
    if(planes == 2 && fPlane > 0){
      this->Pad()->SetTopMargin   (0.110);
      this->Pad()->SetBottomMargin(0.005);
    }
    else if(planes > 2){
      if(fPlane == 1){
	this->Pad()->SetTopMargin   (0.055);
	this->Pad()->SetBottomMargin(0.055);
      }
      else if(fPlane == 2){
	this->Pad()->SetTopMargin   (0.110);
	this->Pad()->SetBottomMargin(0.005);
      }
    }

    TString planeNo = "fTWirePlane";
    planeNo += fPlane;


    TString xtitle = ";Induction Wire;t (tdc)";
    if(geo->Plane(fPlane).SignalType() == geo::kCollection) xtitle = ";Collection Wire;t (tdc)";

    fXLo = -0.005*(geo->Nwires(fPlane)-1);
    fXHi =  1.005*(geo->Nwires(fPlane)-1);
    fYLo = -0.005*this->RawDataDraw()->TotalClockTicks();
    fYHi =  1.01*this->RawDataDraw()->TotalClockTicks();

    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    fOri = rawopt->fAxisOrientation;
    if(fOri > 0){
      fYLo = -0.005*(geo->Nwires(fPlane)-1);
      fYHi =  1.005*(geo->Nwires(fPlane)-1);
      fXLo = -0.005*this->RawDataDraw()->TotalClockTicks();
      fXHi =  1.01*this->RawDataDraw()->TotalClockTicks();
      xtitle = ";t (tdc);InductionWire";
      if(geo->Plane(fPlane).SignalType() == geo::kCollection) xtitle = ";t (tdc);Collection Wire";
    }      

    // make the range of the histogram be the biggest extent 
    // in both directions and then use SetRangeUser() to shrink it down
    // that will allow us to change the axes on the fly
    double min = std::min(fXLo, fYLo);
    double max = std::max(fXHi, fYHi);

    fHisto = new TH1F(*(fPad->DrawFrame(min, min, max, max)));

    fHisto->SetTitleOffset(0.5,"Y");
    fHisto->SetTitleOffset(0.75,"X");
    fHisto->GetYaxis()->SetRangeUser(fYLo, fYHi);
    fHisto->GetYaxis()->SetLabelSize(0.05);
    fHisto->GetYaxis()->CenterTitle();
    fHisto->GetXaxis()->SetRangeUser(fXLo, fXHi);
    fHisto->GetXaxis()->SetLabelSize(0.05);
    fHisto->GetXaxis()->CenterTitle();
    fHisto->Draw("AB");

    fView = new evdb::View2D();
  }

  //......................................................................

  TWireProjPad::~TWireProjPad() 
  {
    if (fHisto) { delete fHisto; fHisto = 0; }
    if (fView)  { delete fView;  fView  = 0; }
  }

  //......................................................................

  void TWireProjPad::Draw(const char* opt) 
  {
    fView->Clear();

    // grab the singleton holding the art::Event
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    if(evt){
      this->SimulationDraw()->MCTruthVectors2D(*evt, fView, fPlane);
      this->RawDataDraw()->   RawDigit2D      (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Wire2D          (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Hit2D           (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Cluster2D       (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Prong2D         (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Vertex2D        (*evt, fView, fPlane);
      this->RecoBaseDraw()->  Event2D         (*evt, fView, fPlane);
    } // if (evt)

    // art::ServiceHandle<evd::GeometryDrawingOptions> drawopt;
    fPad->Clear();
    fPad->cd();

       
    // check if we need to swap the axis ranges
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    if(fOri != rawopt->fAxisOrientation){
      fOri = rawopt->fAxisOrientation;
      double max = fXHi;
      double min = fXLo;
      fXHi = fYHi;
      fXLo = fYLo;
      fYHi = max;
      fYLo = min;
     
      fHisto->GetXaxis()->SetRangeUser(fXLo,fXHi);
      fHisto->GetYaxis()->SetRangeUser(fYLo,fYHi);

      TString xtitle = fHisto->GetXaxis()->GetTitle();
      fHisto->GetXaxis()->SetTitle(fHisto->GetYaxis()->GetTitle());
      fHisto->GetYaxis()->SetTitle(xtitle);
    }

    if (fPlane > 0) fHisto->Draw("X+");
    else            fHisto->Draw("");

      
    // Check if we should zoom the displays
    if (opt==0) {
//       if (drawopt->fAutoZoom) this->AutoZoom();
//       else                    this->ShowFull();
      this->ShowFull();
    }
    
    fPad->GetListOfPrimitives()->FindObject("hframe")->SetBit(TBox::kCannotMove,true);
    fView->Draw();
  }

  //......................................................................

  ///
  /// Automatically zoom the view to a size just larger than the
  /// events. Also ensures that the aspect ratio is the same for the XZ
  /// and YZ projections.
  ///
//   void TWireProjPad::AutoZoom()
//   {
//     double xmin, ymin, zmin;
//     double xmax, ymax, zmax;
//     this->RawDataDraw()->GetLimits(&xmin, &xmax, 
// 				   &ymin, &ymax, 
// 				   &zmin, &zmax);
//     double dx = xmax-xmin;
//     double dy = ymax-ymin;
//     double dz = zmax-zmin;
  
//     if (dx<dy) dx = dy;
//     else       dy = dx;
//     xmin = 0.5*(xmin+xmax)-0.6*dx;
//     xmax = 0.5*(xmin+xmax)+0.6*dx;
//     ymin = 0.5*(ymin+ymax)-0.6*dy;
//     ymax = 0.5*(ymin+ymax)+0.6*dy;
//     zmin -= 0.1*dz;
//     zmax += 0.1*dz;
  
//     fHisto->GetXaxis()->SetRangeUser(zmin,zmax);
//     if (fXorY==kX) fHisto->GetYaxis()->SetRangeUser(xmin,xmax);
//     else           fHisto->GetYaxis()->SetRangeUser(ymin,ymax);
//   }

  //......................................................................
// the override parameter is needed to unzoom to full range when the fAutoZoomInterest is on. 

  void TWireProjPad::ShowFull(int override)
  {
    art::ServiceHandle<geo::Geometry> g;
    
    // x values are wire numbers, y values are ticks of the clock
    int xmin = fXLo;
    int xmax = fXHi;
    int ymax = fYHi;
    int ymin = fYLo;
    
    art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    
    
    if(evdlayoutopt->fAutoZoomInterest && !override)
      {int test=0;
      if(rawopt->fDrawRawDataOrCalibWires == 0)
       test=RawDataDraw()->GetRegionOfInterest((int)fPlane,xmin,xmax,ymin,ymax);
      else
       test=RecoBaseDraw()->GetRegionOfInterest((int)fPlane,xmin,xmax,ymin,ymax);
    
       if(test==-1) return;     
      }
      
    fHisto->GetXaxis()->SetRangeUser(xmin,xmax);
    fHisto->GetYaxis()->SetRangeUser(ymin,ymax);
   
  }

  //......................................................................

  void TWireProjPad::GetWireRange(int* i1, int* i2) const 
  {
    if(fOri < 1){
      *i1 = fHisto->GetXaxis()->GetFirst();
      *i2 = fHisto->GetXaxis()->GetLast();
    }
    else{
      *i1 = fHisto->GetYaxis()->GetFirst();
      *i2 = fHisto->GetYaxis()->GetLast();
    }

  }

  //......................................................................

  void TWireProjPad::SetWireRange(int i1, int i2)
  {
    if(fOri < 1)
    {fHisto->GetXaxis()->SetRange(i1,i2);
     }
    else
    {fHisto->GetYaxis()->SetRange(i1,i2);
    }
  }



  void TWireProjPad::SetZoomRange(int i1, int i2,int y1, int y2)
  {
    
    fHisto->GetXaxis()->SetRangeUser(i1,i2);
    fHisto->GetYaxis()->SetRangeUser(y1,y2);
  }



}// namespace
////////////////////////////////////////////////////////////////////////





