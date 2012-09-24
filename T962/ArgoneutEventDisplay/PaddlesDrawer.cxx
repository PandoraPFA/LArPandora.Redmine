#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TPolyMarker3D.h"

#include "Geometry/Geometry.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/View3D.h"
#include "T962/ArgoneutEventDisplay/ArgoneutDrawingOptions.h"
#include "T962/ArgoneutEventDisplay/PaddlesDrawer.h"
#include "EventDisplay/RecoDrawingOptions.h"

#include "T962/T962_Objects/Paddles.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace argoevd{

  //......................................................................
  PaddlesDrawer::PaddlesDrawer()
  {
  }

  //......................................................................
  PaddlesDrawer::~PaddlesDrawer()
  {
  }

  //......................................................................
  void PaddlesDrawer::DrawPaddlesInfo3D(const art::Event& evt,
                                        evdb::View3D*     view)
  {
    art::ServiceHandle<argoevd::ArgoneutDrawingOptions> argoopt;
    
    art::Handle< t962::Paddles > paddleshandle;
    evt.getByLabel(argoopt->fPaddlesLabel, paddleshandle);
    if(paddleshandle.failedToGet()){
      mf::LogWarning("GetPaddles") << "No Paddles information found.  Skipping.\n";
      return;
    }
    
    int color[4] = {kGray, kGray, kGray, kGray};
    this->FindCoincidences(paddleshandle, color, argoopt->fCoincidenceTime);
    
    //this->GetPaddles(evt, argoopt->fMINOSLabel, paddles); 
    this->PaddlesOutline3D(view, 1, color[0]);
    this->PaddlesOutline3D(view, 2, color[1]);
    this->PaddlesOutline3D(view, 3, color[2]);
    this->PaddlesOutline3D(view, 4, color[3]);
    
  }


   
   //......................................................................
  void PaddlesDrawer::PaddlesOutline3D(evdb::View3D*        view,
                                       int paddlenumber,
                                       int color)
  {
    art::ServiceHandle<geo::Geometry> geo;

    //paddle coordinates
    double length = 58.4;

    double width_offset  = (length - 2.*geo->DetHalfWidth())/2.0;
    double height_offset = (length - 2.*geo->DetHalfHeight())/2.0;
    double length_offset = 30.0;

    
    //tpc active coordinates
    double xlo =  0.                     - width_offset;
    double xhi =  2.*geo->DetHalfWidth() + width_offset;
    double ylo = -geo->DetHalfHeight()   - height_offset;
    double yhi =  geo->DetHalfHeight()   + height_offset;
    double zlo, zhi;
    if(paddlenumber==1){
      zlo =  0.0 - length_offset;
      zhi =  zlo+2.0;
    }
    if(paddlenumber==2){
      zlo =  0.0 - length_offset+8.0;
      zhi =  zlo+2.0;
    }
    if(paddlenumber==3){
      zlo =  geo->DetLength() + length_offset-10.0;
      zhi =  zlo+2.0;
    }
    if(paddlenumber==4){
      zlo =  geo->DetLength() + length_offset-2.0;
      zhi =  zlo+2.0;
    }

    
    int c = color;
    int s = 1;
    int w = 1;
    TPolyLine3D& top = view->AddPolyLine3D(5, c, w, s);
    top.SetPoint(0, xlo, yhi, zlo);
    top.SetPoint(1, xhi, yhi, zlo);
    top.SetPoint(2, xhi, yhi, zhi);
    top.SetPoint(3, xlo, yhi, zhi);
    top.SetPoint(4, xlo, yhi, zlo);

    TPolyLine3D& side = view->AddPolyLine3D(5, c, w, s);
    side.SetPoint(0, xhi, yhi, zlo);
    side.SetPoint(1, xhi, ylo, zlo);
    side.SetPoint(2, xhi, ylo, zhi);
    side.SetPoint(3, xhi, yhi, zhi);
    side.SetPoint(4, xhi, yhi, zlo);

    TPolyLine3D& side2 = view->AddPolyLine3D(5, c, w, s);
    side2.SetPoint(0, xlo, yhi, zlo);
    side2.SetPoint(1, xlo, ylo, zlo);
    side2.SetPoint(2, xlo, ylo, zhi);
    side2.SetPoint(3, xlo, yhi, zhi);
    side2.SetPoint(4, xlo, yhi, zlo);

    TPolyLine3D& bottom = view->AddPolyLine3D(5, c, w, s);
    bottom.SetPoint(0, xlo, ylo, zlo);
    bottom.SetPoint(1, xhi, ylo, zlo);
    bottom.SetPoint(2, xhi, ylo, zhi);
    bottom.SetPoint(3, xlo, ylo, zhi);
    bottom.SetPoint(4, xlo, ylo, zlo);
  
  }

  //......................................................................
  void PaddlesDrawer::FindCoincidences(const art::Handle < t962::Paddles> p,
                                       int* color,
                                       int tol)
  {
    std::cout << *p << std::endl;
    int tdc[4];
    bool val[4] = {false, false, false, false};//bools for Hits within NuMI window.
    
    for(int i = 0; i<4; ++i){
      tdc[i] = p->gettdc(i,0);
      if(tdc[i]<0) tdc[i] = 0;//negative values are garbage
      //Only times of interest are those between 214000-225000,
      //which is approximately the NuMI spill window.
      //Zero out anything outside of this window.
      if(tdc[i]<214000 || tdc[i]>225000) tdc[i] = 0;
      if(tdc[i]!=0) val[i] = true;
    }

    if(!val[0] && !val[1] && !val[2]  && !val[3]) return;//no valid hits

    int d12 = fabs(tdc[0] - tdc[1]);
    int d13 = fabs(tdc[0] - tdc[2]);
    int d14 = fabs(tdc[0] - tdc[3]);
    int d23 = fabs(tdc[1] - tdc[2]);
    int d24 = fabs(tdc[1] - tdc[3]);
    int d34 = fabs(tdc[2] - tdc[3]);
    
    if(val[0]){//paddle1 has a valid hit
      color[0] = kOrange+7;
      if(val[1] && d12<tol) color[1] = kOrange+7;
      if(val[2] && d13<tol) color[2] = kOrange+7;
      if(val[3] && d14<tol) color[3] = kOrange+7;
    }
    if(val[1]){//paddle2 has a valid hit
      if(color[1]==kGray) color[1] = kYellow;
      if(val[2] && d23<tol && color[2]==kGray) color[2] = kYellow;
      if(val[3] && d24<tol && color[3]==kGray) color[3] = kYellow;
    }
    if(val[2]){//paddle3 has a valid hit
      if(color[2]==kGray) color[2] = kGreen;
      if(val[3] && d34<tol && color[3]==kGray) color[3]=kGreen;
    }
    if(val[3]){//paddle4 has a valid hit
      if(color[3]==kGray) color[3]=kMagenta;
    }
  
    return;
    
    
  }

  

}// namespace
////////////////////////////////////////////////////////////////////////
