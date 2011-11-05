/// \file    geo.h
/// \brief   Collect all the geometry header files together
/// \author  brebel@fnal.gov
/// \version $Id: geo.h,v 1.3 2009/04/03 18:13:09 brebel Exp $
#ifndef GEO_GEO_H
#define GEO_GEO_H
#include "Geometry/Geometry.h"
#include "Geometry/WireGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/TPCGeo.h"

/// Detector geometry definition and interface
namespace geo {
  void ProjectToBoxEdge(const double xyz[],
			const double dxyz[],
			double xlo, double xhi,
			double ylo, double yhi,
			double zlo, double zhi,
			double xyzout[]);
  double ClosestApproach(const double point[],
			 const double intercept[],
			 const double slopes[],
			 double closest[]);
  bool CrossesBoundary ( double x0[],
                         double gradient[],
                         double x_lo, double x_hi,
                         double y_lo, double y_hi,
                         double z_lo, double z_hi,
                         double point[] );
}
#endif
