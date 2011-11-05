////////////////////////////////////////////////////////////////////////
/// \file PlaneGeo.cxx
///
/// \version $Id: PlaneGeo.cxx,v 1.12 2010/03/05 19:47:51 bpage Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>


// ROOT includes
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// LArSoft includes
#include "Geometry/geo.h"

namespace geo{

  bool sortByZPos(WireGeo* w1, WireGeo* w2){
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};

    w1->GetCenter(xyz1); w2->GetCenter(xyz2);

    return xyz1[2] < xyz2[2]; 
  }
  //......................................................................

  PlaneGeo::PlaneGeo(std::vector<const TGeoNode*>& path, int depth)
  {
    // build a matrix to take us from the local to the world coordinates
    // in one step
    fGeoMatrix = new TGeoHMatrix(*path[0]->GetMatrix());
    for(int i = 1; i <= depth; ++i){
      fGeoMatrix->Multiply(path[i]->GetMatrix());
    }
  
    // find the wires for the plane so that you can use them later
    this->FindWire(path, depth);

    // sort the wires according to z position
    std::sort(fWire.begin(), fWire.end(), sortByZPos);
  
    // Get the angle of a wire in the plane WRT master z axis
    double theta = fWire[0]->ThetaZ();

    LOG_DEBUG("PlaneGeo") <<"theta = "<< theta << std::endl;
  
    //This uniquely assigns a view to each of the planes
    //in Argoneut, Bo and Microboone, but will need to be made more
    //generic as more geometries are added.
    if(theta >=0.25*TMath::Pi() &&theta <= 0.583333*TMath::Pi()) fView = kU;
    else if(theta > 0.583333*TMath::Pi())                        fView = kV;
    else if(theta < 0.25*TMath::Pi())                            fView = kW;
    else throw cet::exception("BadPlaneView") << "Could not determine view for plane";

    // set all planes to be kVertical.  Bo Geometry is being remade so that 
    // this will always be true. 8/3/11 bjr
    fOrientation = kVertical;

  }

  //......................................................................

  PlaneGeo::~PlaneGeo()
  {
    for(unsigned int i = 0; i < fWire.size(); ++i)
      if(fWire[i]) delete fWire[i];
  
    fWire.clear();

    if(fGeoMatrix) delete fGeoMatrix;

  }

  //......................................................................

  void PlaneGeo::FindWire(std::vector<const TGeoNode*>& path,
			  unsigned int depth) 
  {
    // Check if the current node is a wire
    const char* wire = "volTPCWire";
    if(strncmp(path[depth]->GetName(), wire, strlen(wire)) == 0){
      this->MakeWire(path, depth);
      return;
    }
  
    // Explore the next layer down
    unsigned int deeper = depth+1;
    if (deeper>=path.size()) {
      throw cet::exception("ExceededMaxDepth") << "Exceeded maximum depth";
    }
    const TGeoVolume* v = path[depth]->GetVolume();
    int nd = v->GetNdaughters();
    for (int i=0; i<nd; ++i) {
      path[deeper] = v->GetNode(i);
      this->FindWire(path, deeper);
    }
  }

  //......................................................................

  void PlaneGeo::MakeWire(std::vector<const TGeoNode*>& path, int depth) 
  {
    fWire.push_back(new WireGeo(path, depth));
  }

  //......................................................................

  void PlaneGeo::LocalToWorld(const double* plane, double* world) const
  {
    fGeoMatrix->LocalToMaster(plane, world);
  }

  //......................................................................

  void PlaneGeo::LocalToWorldVect(const double* plane, double* world) const
  {
    fGeoMatrix->LocalToMasterVect(plane, world);
  }

  //......................................................................

  void PlaneGeo::WorldToLocal(const double* world, double* plane) const
  {
    fGeoMatrix->MasterToLocal(world, plane);
  }

  //......................................................................

  const TVector3 PlaneGeo::WorldToLocal( const TVector3& world ) const
  {
    double worldArray[4];
    double localArray[4];
    worldArray[0] = world.X();
    worldArray[1] = world.Y();
    worldArray[2] = world.Z();
    worldArray[3] = 1.; 
    fGeoMatrix->MasterToLocal(worldArray,localArray);
    return TVector3(localArray);
  }

  //......................................................................

  const TVector3 PlaneGeo::LocalToWorld( const TVector3& local ) const
  {
    double worldArray[4];
    double localArray[4];
    localArray[0] = local.X();
    localArray[1] = local.Y();
    localArray[2] = local.Z();
    localArray[3] = 1.;
    fGeoMatrix->LocalToMaster(localArray,worldArray);
    return TVector3(worldArray);
  }

  //......................................................................

  // Convert a vector from world frame to the local plane frame
  // \param world : 3-D array. Vector in world coordinates; input.
  // \param plane : 3-D array. Vector in plane coordinates; plane.
  void PlaneGeo::WorldToLocalVect(const double* world, double* plane) const
  {
    fGeoMatrix->MasterToLocalVect(world,plane);
  }

}
////////////////////////////////////////////////////////////////////////
