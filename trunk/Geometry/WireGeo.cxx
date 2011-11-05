////////////////////////////////////////////////////////////////////////
/// \file  WireGeo.cxx
/// \brief Encapsulate the geometry of a wire
///
/// \version $Id: WireGeo.cxx,v 1.8 2010/03/05 05:30:56 bpage Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "Geometry/geo.h"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TMath.h"

namespace geo{

  /// Construct a wire geometry
  /// \param path  : List of TGeoNodes that take us to this wire
  /// \param depth : Size of "path" list

  //-----------------------------------------
  WireGeo::WireGeo(std::vector<const TGeoNode*>& path, int depth) 
  {
    fWireNode = path[depth];

    /// uncomment the following to check the paths to the wires
    ///   std::string p(base);
    ///   for(int i = 0; i <= depth; ++i){
    ///     p += "/";
    ///     p += path[i]->GetName();
    ///   }
    ///   std::cout << p.c_str() << std::endl;
  
    // Build the matrix that takes us to the top world frame
    // build a matrix to take us from the local to the world coordinates
    // in one step
    fGeoMatrix = new TGeoHMatrix(*path[0]->GetMatrix());
    for(int i = 1; i <= depth; ++i){
      fGeoMatrix->Multiply(path[i]->GetMatrix());
    }
  }

  //......................................................................
  WireGeo::~WireGeo()
  {
    //if(fGeoMatrix) delete fGeoMatrix;
  
    return;
  }  

  //......................................................................

  /// Transform a position from local frame to world frame
  /// \param local : 3D array. Position in the local frame  Input.
  /// \param world : 3D array. Position in the world frame. Returned.
  void WireGeo::LocalToWorld(const double* local, double* world) const
  {
    fGeoMatrix->LocalToMaster(local,world);
  }

  //......................................................................    

  /// Transform a 3-vector from local frame to world frame
  /// \param local : 3D array. Position in the local frame  Input.
  /// \param world : 3D array. Position in the world frame. Returned.
  void WireGeo::LocalToWorldVect(const double* local, double* world) const
  {
    fGeoMatrix->LocalToMasterVect(local,world);
  }
    
  //......................................................................

  /// Transform a position from world frame to local frame
  /// \param world : 3D array. Position in the world frame. Input.
  /// \param local : 3D array. Position in the local frame  Returned.
  void WireGeo::WorldToLocal(const double* local, double* world) const
  {
    fGeoMatrix->MasterToLocal(local,world);
  }

  //......................................................................

  /// Transform a 3-vector from world frame to local frame
  /// \param world : 3D array. Position in the world frame. Input.
  /// \param local : 3D array. Position in the local frame  Returned.
  void WireGeo::WorldToLocalVect(const double* local, double* world) const
  {
    fGeoMatrix->MasterToLocalVect(local,world);
  }

  //......................................................................

  /// Return the center position of a wire.
  /// \param xyz : 3-D array. The returned location.
  /// \param localz : Distance along the length of the wire
  /// (cm). Default is center of wire
  void WireGeo::GetCenter(double* xyz, double localz) const
  {
    double xyzLocal[3] = {0.,0.,localz};
    this->LocalToWorld(xyzLocal, xyz);
  }

  //......................................................................

  double WireGeo::RMax() const 
  {
    return ((TGeoTube*)fWireNode->GetVolume()->GetShape())->GetRmax();
  }
    
  //......................................................................

  double WireGeo::HalfL() const 
  {
    return ((TGeoTube*)fWireNode->GetVolume()->GetShape())->GetDZ();
  }

  //......................................................................

  double WireGeo::RMin() const 
  {
    return ((TGeoTube*)fWireNode->GetVolume()->GetShape())->GetRmin();
  }

  //......................................................................
  double WireGeo::ThetaZ(bool degrees) const
  {
    static double xyzCenter[3] = {0,0,0};
    static double xyzEnd[3] = {0,0,0};
    static double halfL =this->HalfL();
    this->GetCenter(xyzCenter,0.0);
    this->GetCenter(xyzEnd,halfL);
    //either y or x will be 0, so ading both will always catch the right
    //one
    double angle = (xyzEnd[1]-xyzCenter[1]+xyzEnd[0]-xyzCenter[0]) / 
      TMath::Abs(xyzEnd[1]-xyzCenter[1]+xyzCenter[0]-xyzEnd[0]) * 
      TMath::ACos((xyzEnd[2] - xyzCenter[2])/halfL);
    //This ensures we are looking at the angle between 0 and Pi
    //as if the wire runs at one angle it also runs at that angle +-Pi
    if(angle < 0) angle +=TMath::Pi(); 
    if(degrees) angle*=180.0/TMath::Pi();
    return angle;
  }

}
////////////////////////////////////////////////////////////////////////
