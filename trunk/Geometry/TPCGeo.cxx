////////////////////////////////////////////////////////////////////////
/// \file TPCGeo.cxx
///
/// \version $Id: TPCGeo.cxx,v 1.12 2010/03/05 19:47:51 bpage Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>


// ROOT includes
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include <TGeoBBox.h>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// LArSoft includes
#include "Geometry/geo.h"

namespace geo{

  //......................................................................
  // Define sort order for detector planes.
  static bool plane_sort(const PlaneGeo* p1, const PlaneGeo* p2) 
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.};
    p1->LocalToWorld(local, xyz1);
    p2->LocalToWorld(local, xyz2);

    // std::cout << xyz1[0] << " " << xyz2[0] << " " << p1->Orientation() << std::endl;
    
    // All planes for all geometries are assumed to be in the yz plane, ie vertical
    // Bo geometry will shortly comply with this assumption. 8/3/11 bjr
    return xyz1[0]>xyz2[0];
  }

  //......................................................................
  TPCGeo::TPCGeo(std::vector<const TGeoNode*>& path, int depth)
  {
    // all planes are going to be contained in the volume named volTPC
    // now get the total volume of the TPC
    TGeoVolume *vc = gGeoManager->FindVolumeFast("volTPC");
    if(vc){
      fTotalVolume = gGeoManager->FindVolumeFast("volTPC");
      if(!vc){ 
	mf::LogWarning("Geometry") << "cannot find detector outline volume - bail ungracefully";
	assert(0);
      }

      TGeoVolume *vca = gGeoManager->FindVolumeFast("volTPCActive");
      if(vca) fActiveVolume = vca;
      else    fActiveVolume = vc;
    }
    mf::LogInfo("Geometry") << "detector total  volume is " << fTotalVolume->GetName();
    mf::LogInfo("Geometry") << "detector active volume is " << fActiveVolume->GetName();

    // build a matrix to take us from the local to the world coordinates
    // in one step
    fGeoMatrix = new TGeoHMatrix(*path[0]->GetMatrix());
    for(int i = 1; i <= depth; ++i){
      fGeoMatrix->Multiply(path[i]->GetMatrix());
    }
  
    // find the wires for the plane so that you can use them later
    this->FindPlane(path, depth);

    // sort the wires according to z position
    std::sort(fPlanes.begin(), fPlanes.end(), plane_sort);

    // set the plane pitch for this TPC
    double origin[3] = {0.};
    double xyz[3]  = {0.};
    fPlanes[0]->LocalToWorld(origin,xyz);
    double xyz1[3] = {0.};
    fPlaneLocation.clear();
    fPlaneLocation.resize(fPlanes.size());
    for(unsigned int i = 0; i < fPlaneLocation.size(); ++i) fPlaneLocation[i].resize(3);
    fPlane0Pitch.clear();
    fPlane0Pitch.resize(this->Nplanes(), 0.);
    for(size_t p = 0; p < this->Nplanes(); ++p){
      fPlanes[p]->SetSignalType(geo::kInduction); //<set all planes to be induction for now
      fPlanes[p]->LocalToWorld(origin,xyz1);
      if(p > 0) fPlane0Pitch[p] = fPlane0Pitch[p-1] + fabs(xyz1[0]-xyz[0]);
      else      fPlane0Pitch[p] = 0.;
      xyz[0] = xyz1[0];
      fPlaneLocation[p][0] = xyz1[0];
      fPlaneLocation[p][1] = xyz1[1];
      fPlaneLocation[p][2] = xyz1[2];
    }

    // now set the collection plane type
    fPlanes[fPlanes.size()-1]->SetSignalType(geo::kCollection); 

    unsigned int chan = 0;
    for(unsigned int p = 0; p < fPlanes.size(); ++p){
      for(unsigned int w = 0; w < fPlanes[p]->Nwires(); ++w){
	PlaneWirePair pwp(p, w);
	fChannelMap[chan] = pwp;
	++chan;
      }
    }

    // determine the drift direction of the electrons in the TPC
    // first get the location of the planes in the world coordinates
    double planeworld[3] = {0.};
    double tpcworld[3] = {0.};
    
    fPlanes[0]->LocalToWorld(origin, planeworld);

    // now get the origin of the TPC in world coordinates
    this->LocalToWorld(origin, tpcworld);

    // check to see if the x coordinates change between the tpc
    // origin and the plane origin, and if so in which direction
    fDriftDirection = geo::kUnknownDrift;
    if     ( tpcworld[0] > 1.01*planeworld[0] ) fDriftDirection = geo::kNegX;
    else if( tpcworld[0] < 0.99*planeworld[0] ) fDriftDirection = geo::kPosX;
    else if( tpcworld[1] < 0.99*planeworld[1] ) fDriftDirection = geo::kPosY;

    return;
  }

  //......................................................................
  TPCGeo::~TPCGeo()
  {
    for(unsigned int i = 0; i < fPlanes.size(); ++i)
      if(fPlanes[i]) delete fPlanes[i];
  
    fPlanes.clear();

    if(fGeoMatrix)    delete fGeoMatrix;
    if(fActiveVolume) delete fActiveVolume;
  }

  //......................................................................
  void TPCGeo::FindPlane(std::vector<const TGeoNode*>& path,
			 unsigned int depth) 
  {

    const char* nm = path[depth]->GetName();
    if( (strncmp(nm, "volTPCPlane", 11) == 0) ){
      this->MakePlane(path,depth);
      return;
    }

    //explore the next layer down
    unsigned int deeper = depth+1;
    if(deeper >= path.size()){
      throw cet::exception("BadTGeoNode") << "exceeded maximum TGeoNode depth";
    }

    const TGeoVolume *v = path[depth]->GetVolume();
    int nd = v->GetNdaughters();
    for(int i = 0; i < nd; ++i){
      path[deeper] = v->GetNode(i);
      this->FindPlane(path, deeper);
    }
  
  }

  //......................................................................
  void TPCGeo::MakePlane(std::vector<const TGeoNode*>& path, int depth) 
  {
    fPlanes.push_back(new PlaneGeo(path, depth));
  }

  //......................................................................
  const PlaneGeo& TPCGeo::Plane(unsigned int iplane) const
  {
    if(iplane >= fPlanes.size()){
      mf::LogWarning("PlaneOutOfRange") << "Request for non-existant plane " << iplane
					<< " bail";
      assert(0);
    }

    return *fPlanes[iplane];
  }

  //......................................................................
  double TPCGeo::ActiveHalfWidth()  const 
  {
    return ((TGeoBBox*)fActiveVolume->GetShape())->GetDX();
  }

  //......................................................................
  double TPCGeo::ActiveHalfHeight() const 
  {
    return ((TGeoBBox*)fActiveVolume->GetShape())->GetDY();
  }

  //......................................................................
  double TPCGeo::ActiveLength() const
  { 
    return 2.0*((TGeoBBox*)fActiveVolume->GetShape())->GetDZ();
  }

  //......................................................................
  double TPCGeo::HalfWidth()  const 
  {
    return ((TGeoBBox*)fTotalVolume->GetShape())->GetDX();
  }

  //......................................................................
  double TPCGeo::HalfHeight() const 
  {
    return ((TGeoBBox*)fTotalVolume->GetShape())->GetDY();
  }

  //......................................................................
  double TPCGeo::Length() const
  { 
    return 2.0*((TGeoBBox*)fTotalVolume->GetShape())->GetDZ();
  }

  //......................................................................
  // returns distance between plane 0 to each of the remaining planes 
  // not the distance between two consecutive planes  
  double TPCGeo::Plane0Pitch(unsigned int p) const
  {
    return fPlane0Pitch[p];
  }

  //......................................................................
  // returns xyz location of planes in TPC
  const double* TPCGeo::PlaneLocation(unsigned int p) const
  {
    return &fPlaneLocation[p][0];
  }

  //......................................................................
  double TPCGeo::PlanePitch(unsigned int p1, 
                            unsigned int p2) const
  {
    double xyz[3] = {0.};
    double xyz1[3] = {0.};

    //every plane has a wire 0
    this->Plane(p1).Wire(0).GetCenter(xyz);
    this->Plane(p2).Wire(0).GetCenter(xyz1);
    
    return fabs(xyz1[0]-xyz[0]);
  }

  //......................................................................
  // This method returns the distance between the specified wires.
  // w1 < w2.  The wires are assumed to be on the same plane
  double TPCGeo::WirePitch(unsigned int w1,  
			   unsigned int w2,  
			   unsigned int plane) const
  { 
    double xyz[3] = {0.};
    double xyz1[3] = {0.};

    this->Plane(plane).Wire(w1).GetCenter(xyz,  this->Plane(plane).Wire(w1).HalfL());
    this->Plane(plane).Wire(w2).GetCenter(xyz1, this->Plane(plane).Wire(w2).HalfL());

    if(xyz1[2] - xyz[2] < 0.01){
      this->Plane(plane).Wire(w1).GetCenter(xyz,  -this->Plane(plane).Wire(w1).HalfL());
      this->Plane(plane).Wire(w2).GetCenter(xyz1, -this->Plane(plane).Wire(w2).HalfL());
    }

    double thetaz = this->Plane(plane).Wire(w2).ThetaZ();
    double pitch = fabs((xyz1[2]-xyz[2])*TMath::Sin(thetaz) 
			-(xyz1[1]-xyz[1])*TMath::Cos(thetaz));
    
    return pitch;
  }

  //......................................................................
  unsigned int TPCGeo::PlaneWireToTPCChannel(unsigned int plane,
					     unsigned int wire) const
  {
    PlaneWirePair pwp(plane,wire);
    TPCChannelMap::const_iterator itr = fChannelMap.begin();
    while(itr != fChannelMap.end() ){
      if(itr->second == pwp) return itr->first;
      itr++;
    }

    // if we got here the plane and wire pair are not in the tpc
    throw cet::exception("TPCGeo::PlaneWireToTPCChannel: NO CHANNEL FOUND ") << "for plane,wire: "
									     << plane << "," << wire
									     << "\n returning UINT_MAX";
    return UINT_MAX;
  }

  //......................................................................
  const WireGeo& TPCGeo::ChannelToWire(unsigned int  channel,
				       unsigned int &plane,
				       unsigned int &wire) const
  {
    //the conversion from channel number to plane and wire is straightforward
    TPCChannelMap::const_iterator itr = fChannelMap.find(channel);
    if( itr == fChannelMap.end() ){
      throw cet::exception("geo::ChannelToWire(): BAD CHANNEL Number ") << "channel " << channel << " "
									<< " not found in map";
      
    }

    plane = itr->second.first;
    wire  = itr->second.second;

    return this->Plane(plane).Wire(wire);
  }

  //......................................................................
  unsigned int TPCGeo::NearestChannel(double* worldPos) const
  {
    // This routine is the most frequently called in the entire
    // simulation, during the voxel->electron cluster calculation. Try
    // to make this code as efficient as possible.

    // First, a bit of "buffering": This routine is usually called in
    // succession for different electron clusters generated from the
    // same voxel. This means that the value of worldPos may not
    // change much, and therefore neither would the result. So test if
    // worldPos for this call of the routine is close enough to the
    // last one for which we did a full calculation; if it is, skip
    // the calculation and return the previous result.
    
    // Get every little bit of speed by not using an array for the
    // previous point.
    static double lastPos0, lastPos1, lastPos2;
    static unsigned int nearest = 0;
    static double closeEnough = 0.;
    static bool firstCalculation = true;

    if ( firstCalculation ){
      firstCalculation = false;
      // What do we mean by "close enough" to the previous point?
      // For now, let's take it to be one-tenth the of the wire
      // spacing; this is typically the same size as the voxels used
      // in the simulation.
      closeEnough = this->WirePitch( 0, 1, 0 ) / 10.;
    }
    else{
      if ( fabs(worldPos[0] - lastPos0) < closeEnough  &&
	   fabs(worldPos[1] - lastPos1) < closeEnough  &&
	   fabs(worldPos[2] - lastPos2) < closeEnough  )
	{ return nearest; }
    }

    // If we get to this line, we're doing the full calculation. Save
    // the current point.
    lastPos0 = worldPos[0];
    lastPos1 = worldPos[1];
    lastPos2 = worldPos[2];

    // To speed up the calculation, pre-compute as much constant 
    // information as we can. 

    // Number of planes.
    static int nplanes;
    // x-coordinate of each plane.
    static std::vector<double> planex;
    // Number of wires in each plane.
    static std::vector<int> nwires;
    // For each plane, for each wire, we need to store some information.
    // For wireLow and wireVector, we're storing three-vectors.
    static std::vector< std::vector< std::vector< double > > > wireLow;
    static std::vector< std::vector< std::vector< double > > > wireVect;
    static std::vector< std::vector< double > > wireMag;
    static std::vector< std::vector< unsigned int > > planeWireToChannel;

    static bool initialized = false;
    if ( ! initialized ){
      initialized = true;
      
      // The zero three-vector:
      const std::vector< double > zeroVector(3,0.);
      
      // std::cout<<"number of planes is  "<<this->Nplanes()<<std::endl;
      nplanes = this->Nplanes();
      // Now that we know the number of planes, allocate storage.
      planex.resize(nplanes);
      nwires.resize(nplanes);
      wireLow.resize(nplanes);
      wireVect.resize(nplanes);
      wireMag.resize(nplanes);
      planeWireToChannel.resize(nplanes);
      
      // Dummy arrays for coordinate transformation.
      double pos[3] = {0.};
      static const double origin[3] = {0.};
      for(int p = 0; p < nplanes; ++p){
	
	// Calculate the x-coordinate of the plane in the world
	// volume.
	this->Plane(p).LocalToWorld(origin, pos);
	planex[p] = pos[0];
	
	// Wire values.
	nwires[p] = this->Plane(p).Nwires();
	// Pre-allocate the length of the array for each plane.
	std::vector< std::vector< double > > lowPoints(nwires[p],zeroVector);
	std::vector< std::vector< double > > wireCenters(nwires[p],zeroVector);
	std::vector< double > mags(nwires[p]);
	std::vector< unsigned int > wireToChannel(nwires[p]);
	
	for ( int w = 0; w < nwires[p]; ++w ){
	  const WireGeo wg = this->Plane(p).Wire(w);
	  
	  // Again, dummy arrays for geometry routines.
	  static double pos1[3] = {0.};
	  static double pos2[3] = {0.};
	  // get two points on the wire in order to form a vector - use the end points
	  wg.GetCenter(pos1, -wg.HalfL());
	  wg.GetCenter(pos2,  wg.HalfL());
	  double mag2 = 0.;
	  for ( unsigned int c=0; c < 3; c++ ){
	    lowPoints[w][c] = pos1[c];
	    wireCenters[w][c] = pos2[c] - lowPoints[w][c];
	    mag2 += wireCenters[w][c] * wireCenters[w][c];
	  }
	  mags[w] = sqrt( mag2 );
	  wireToChannel[w] = this->PlaneWireToTPCChannel(p, w);
	  
	} // loop over wires
	
	// Save the wire information for each plane.
	wireLow[p] = lowPoints;
	wireVect[p] = wireCenters;
	wireMag[p] = mags;
	planeWireToChannel[p] = wireToChannel;
	
      } // end loop over planes
    } // end initialization
    
    
    //the planes increase in number with decreasing x value
    //if x is between two planes, take the larger plane number
    int plane = 0;
    for(int i = 0; i < nplanes - 1; ++i){
      
      // Small offset since testing real numbers for equality is always tricky.
      static const double epsilon = 6e-8;
      
      if((worldPos[0] + epsilon) < planex[i]
	 && (worldPos[0] + epsilon) > planex[i+1]){ 
	plane = i+1;
      }
      
    } //end loop over planes
    

    int wire = 0;
    
    double minDist = 1.e12;
    double dist = 0.;
    double lastDist = 0.;
    // Squeeze every cycle of performance by making these
    // simple variables instead of arrays.
    static double worldPosVect0, worldPosVect1, worldPosVect2;
    static double crossProd0, crossProd1, crossProd2;
    static double mag2;

    while( wire < nwires[plane] ){

      // make vectors from the wire center to the test position and from the center to 
      // the point along its length
      worldPosVect0 = worldPos[0] - wireLow[plane][wire][0];
      worldPosVect1 = worldPos[1] - wireLow[plane][wire][1];
      worldPosVect2 = worldPos[2] - wireLow[plane][wire][2];

      crossProd0 = worldPosVect1 * wireVect[plane][wire][2] - worldPosVect2 * wireVect[plane][wire][1];
      crossProd1 = worldPosVect2 * wireVect[plane][wire][0] - worldPosVect0 * wireVect[plane][wire][2];
      crossProd2 = worldPosVect0 * wireVect[plane][wire][1] - worldPosVect1 * wireVect[plane][wire][0];
      
      //distance from the position to the wire 
      mag2 = crossProd0*crossProd0
	+ crossProd1*crossProd1
	+ crossProd2*crossProd2;
      dist = sqrt(mag2)/wireMag[plane][wire];
      
      if(dist < minDist){
	minDist = dist;
	nearest = planeWireToChannel[plane][wire];
      }

      // If the distance to the nearest wire was going down as we scan
      // through the wires, and now it's going up again, then we're
      // done; we've already found the closest wire. (If I were smarter,
      // I'd find a way to do this using a binary search.)
      if ( wire == 0 ) lastDist = dist;
      if ( dist > lastDist ) break; 
      lastDist = dist;

      ++wire;
    }
    
    return nearest;

  }


  //......................................................................
  void TPCGeo::LocalToWorld(const double* tpc, double* world) const
  {
    fGeoMatrix->LocalToMaster(tpc, world);
  }

  //......................................................................
  void TPCGeo::LocalToWorldVect(const double* tpc, double* world) const
  {
    fGeoMatrix->LocalToMasterVect(tpc, world);
  }

  //......................................................................

  void TPCGeo::WorldToLocal(const double* world, double* tpc) const
  {
    fGeoMatrix->MasterToLocal(world, tpc);
  }

  //......................................................................

  const TVector3 TPCGeo::WorldToLocal( const TVector3& world ) const
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

  const TVector3 TPCGeo::LocalToWorld( const TVector3& local ) const
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
  void TPCGeo::WorldToLocalVect(const double* world, double* plane) const
  {
    fGeoMatrix->MasterToLocalVect(world,plane);
  }

}
////////////////////////////////////////////////////////////////////////
