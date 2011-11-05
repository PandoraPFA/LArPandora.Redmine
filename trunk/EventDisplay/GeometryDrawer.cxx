/// \file    GeometryDrawer.cxx
/// \brief   Class to aid in the rendering of Geometry objects
/// \author  messier@indiana.edu
/// \version $Id: GeometryDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"

#include "Geometry/geo.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/View3D.h"
#include "EventDisplay/GeometryDrawer.h"
//#include "EventDisplay/GeometryDrawingOptions.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd{

  //......................................................................
  GeometryDrawer::GeometryDrawer()
  {
  }

  //......................................................................
  GeometryDrawer::~GeometryDrawer()
  {
  }

  //......................................................................
  void GeometryDrawer::DetOutline3D(evdb::View3D*        view)
  {
    art::ServiceHandle<geo::Geometry> geo;

    double xlo =  0.;
    double xhi =  2.*geo->DetHalfWidth();
    double ylo = -geo->DetHalfHeight();
    double yhi =  geo->DetHalfHeight();
    double zlo =  0.0;
    double zhi =  geo->DetLength();
  
    int c = kGray;
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

    c = kGray+2;
    s = 1;
    w = 1;
    double z = zlo;
    // Grid running along x and y at constant z
    for (;;) {
      TPolyLine3D& gridt = view->AddPolyLine3D(2, c, s, w);
      gridt.SetPoint(0, xlo, ylo, z);
      gridt.SetPoint(1, xhi, ylo, z);

      TPolyLine3D& grids = view->AddPolyLine3D(2, c, s, w);
      grids.SetPoint(0, xhi, ylo, z);
      grids.SetPoint(1, xhi, yhi, z);

      z += 10.0;
      if (z>zhi) break;
    }

    // Grid running along z at constant x
    double x = 0.0;
    for (;;) {
      TPolyLine3D& gridt = view->AddPolyLine3D(2, c, s, w);
      gridt.SetPoint(0, x, ylo, zlo);
      gridt.SetPoint(1, x, ylo, zhi);
      x += 10.0;
      if (x>xhi) break;
    }

    // Grid running along z at constant y
    double y = 0.0;
    for (;;) {
      TPolyLine3D& grids = view->AddPolyLine3D(2, c, s, w);
      grids.SetPoint(0, xhi, y, zlo);
      grids.SetPoint(1, xhi, y, zhi);
      y += 10.0;
      if (y>yhi) break;
    }
    y = -10.0;
    for (;;) {
      TPolyLine3D& grids = view->AddPolyLine3D(2, c, s, w);
      grids.SetPoint(0, xhi, y, zlo);
      grids.SetPoint(1, xhi, y, zhi);
      y -= 10.0;
      if (y<ylo) break;
    }

    // Indicate coordinate system
    double x0 = -0.20;     // Center location of the key
    double y0 =  1.10*ylo; // Center location of the key
    double z0 = -0.10*zhi; // Center location of the key  
    double sz =  0.20*zhi; // Scale size of the key in z direction

    c = kBlue;

    TPolyLine3D& xaxis = view->AddPolyLine3D(2, c, s, w);
    TPolyLine3D& yaxis = view->AddPolyLine3D(2, c, s, w);
    TPolyLine3D& zaxis = view->AddPolyLine3D(2, c, s, w);
    xaxis.SetPoint(0, x0,    y0, z0);
    xaxis.SetPoint(1, sz+x0, y0, z0);

    yaxis.SetPoint(0, x0, y0,     z0);
    yaxis.SetPoint(1, x0, y0+sz,  z0);

    zaxis.SetPoint(0, x0, y0, z0);
    zaxis.SetPoint(1, x0, y0, z0+sz);

    TPolyLine3D& xpoint = view->AddPolyLine3D(3, c, s, w);
    TPolyLine3D& ypoint = view->AddPolyLine3D(3, c, s, w);
    TPolyLine3D& zpoint = view->AddPolyLine3D(3, c, s, w);
  
    xpoint.SetPoint(0, 0.95*sz+x0, y0, z0-0.05*sz);
    xpoint.SetPoint(1, 1.00*sz+x0, y0, z0);
    xpoint.SetPoint(2, 0.95*sz+x0, y0, z0+0.05*sz);

    ypoint.SetPoint(0, x0, 0.95*sz+y0, z0-0.05*sz);
    ypoint.SetPoint(1, x0, 1.00*sz+y0, z0);
    ypoint.SetPoint(2, x0, 0.95*sz+y0, z0+0.05*sz);

    zpoint.SetPoint(0, x0-0.05*sz, y0, 0.95*sz+z0);
    zpoint.SetPoint(1, x0+0.00*sz, y0, 1.00*sz+z0);
    zpoint.SetPoint(2, x0+0.05*sz, y0, 0.95*sz+z0);

    TPolyLine3D& zleg = view->AddPolyLine3D(4, c, s, w);  
    zleg.SetPoint(0,  x0-0.05*sz, y0+0.05*sz, z0+1.05*sz);
    zleg.SetPoint(1,  x0+0.05*sz, y0+0.05*sz, z0+1.05*sz);
    zleg.SetPoint(2,  x0-0.05*sz, y0-0.05*sz, z0+1.05*sz);
    zleg.SetPoint(3,  x0+0.05*sz, y0-0.05*sz, z0+1.05*sz);

    TPolyLine3D& yleg = view->AddPolyLine3D(5, c, s, w);  
    yleg.SetPoint(0,  x0-0.05*sz, y0+1.15*sz, z0);
    yleg.SetPoint(1,  x0+0.00*sz, y0+1.10*sz, z0);
    yleg.SetPoint(2,  x0+0.00*sz, y0+1.05*sz, z0);
    yleg.SetPoint(3,  x0+0.00*sz, y0+1.10*sz, z0);
    yleg.SetPoint(4,  x0+0.05*sz, y0+1.15*sz, z0);

    TPolyLine3D& xleg = view->AddPolyLine3D(7, c, s, w);  
    xleg.SetPoint(0,  x0+1.05*sz, y0+0.05*sz, z0-0.05*sz);
    xleg.SetPoint(1,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(2,  x0+1.05*sz, y0+0.05*sz, z0+0.05*sz);
    xleg.SetPoint(3,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(4,  x0+1.05*sz, y0-0.05*sz, z0-0.05*sz);
    xleg.SetPoint(5,  x0+1.05*sz, y0+0.00*sz, z0-0.00*sz);
    xleg.SetPoint(6,  x0+1.05*sz, y0-0.05*sz, z0+0.05*sz);

  }


}// namespace
////////////////////////////////////////////////////////////////////////
