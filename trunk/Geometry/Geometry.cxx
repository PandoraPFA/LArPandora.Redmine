////////////////////////////////////////////////////////////////////////
/// \file Geometry.cxx
///
/// \version $Id: Geometry.cxx,v 1.19 2010/04/27 14:20:10 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// C/C++ includes
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib> // for std::abort()

// lar includes
#include "Geometry/geo.h"
#include "SummaryData/summary.h"

// ROOT includes
#include <TMath.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TVector3.h>
#include <Rtypes.h>

// Framework includes
#include "cetlib/exception.h"
#include "cetlib/search_path.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {

  //......................................................................
  // Define sort order for detector planes.
  static bool tpc_sort(const TPCGeo* t1, const TPCGeo* t2) 
  {
    double xyz1[3] = {0.}, xyz2[3] = {0.};
    double local[3] = {0.};
    t1->LocalToWorld(local, xyz1);
    t2->LocalToWorld(local, xyz2);

    // std::cout << xyz1[0] << " " << xyz2[0] << " " << p1->Orientation() << std::endl;

    return xyz1[0]>xyz2[0];
  }

  //......................................................................
  // Constructor.
  Geometry::Geometry(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
    : fSurfaceY        (pset.get< double          >("SurfaceY"               ))
    , fDetectorName    (pset.get< std::string     >("Name"                   ))
    , fMinWireZDist    (pset.get< double          >("MinWireZDist",     3.0  ))
    , fDisableWiresInG4(pset.get< bool            >("DisableWiresInG4", false))
    , fForceUseFCLOnly (pset.get< bool            >("ForceUseFCLOnly" , false))
    , fNBadVol(0)
  {
    reg.watchPreBeginRun(this, &Geometry::preBeginRun);

    fDetectorName.ToLower();
    if     (fDetectorName.Contains("argoneut")  ) fDetId = geo::kArgoNeuT;
    else if(fDetectorName.Contains("microboone")) fDetId = geo::kMicroBooNE;
    else if(fDetectorName.Contains("lbne")      ) fDetId = geo::kLBNE;
    else if(fDetectorName.Contains("bo")        ) fDetId = geo::kBo;
    
    // Search all reasonable locations for the GDML file that contains
    // the detector geometry.
    // constructor decides if initialized value is a path or an environment variable
    cet::search_path sp("FW_SEARCH_PATH");

    std::string GDMLFileName = pset.get< std::string >("GDML");
    std::string ROOTFileName = pset.get< std::string >("ROOT");

    if(fDisableWiresInG4) GDMLFileName.insert(GDMLFileName.find(".gdml"),"_nowires")  ;
    
    sp.find_file(GDMLFileName, fGDMLfile);
    sp.find_file(ROOTFileName, fROOTfile);

    this->LoadGeometryFile(fGDMLfile, fROOTfile);
  }

  //......................................................................
  Geometry::~Geometry() 
  {
    for (unsigned int i=0; i<fTPCs.size(); ++i) {
      if (fTPCs[i]) { delete fTPCs[i]; fTPCs[i] = 0; }
    }
    
  }

  //......................................................................
  void Geometry::preBeginRun(art::Run const& run)
  {
    // check here to see if we need to load a new geometry.
    // get the detector id from the run object
    std::vector< art::Handle<sumdata::RunData> > rdcol;
    run.getManyByType(rdcol);
    if(rdcol.size() > 0 && !fForceUseFCLOnly){ 

      if(rdcol[0]->DetId() == fDetId) return;
      
      fDetId = rdcol[0]->DetId();

      std::string relpathgdml("Geometry/gdml/");
      std::string relpathroot("Geometry/gdml/");
      
      switch(fDetId){
      case geo::kBo         : 
	relpathgdml += "bo";         relpathroot += "bo.root";         break;
      case geo::kArgoNeuT   : 
	relpathgdml += "argoneut";   relpathroot += "argoneut.root";   break;
      case geo::kMicroBooNE : 
	relpathgdml += "microboone"; relpathroot += "microboone.root"; break;
      case geo::kLBNE       : 
	relpathgdml += "lbne";       relpathroot += "lbne.root";       break;
      default                        : 
	mf::LogWarning("LoadNewGeometry") << "detid invalid, " << fDetId << " give up"; assert(0);
      }
      
      if(fDisableWiresInG4) relpathgdml+="_nowires.gdml";
      else                  relpathgdml+=".gdml";

      
      // constructor decides if initialized value is a path or an environment variable
      cet::search_path sp("FW_SEARCH_PATH");
      
      sp.find_file(relpathgdml, fGDMLfile);
      sp.find_file(relpathroot, fROOTfile);
      
      this->LoadGeometryFile(fGDMLfile, fROOTfile);
    }
    else 
      mf::LogWarning("LoadNewGeometry") << "cannot find sumdata::RunData object to grab detid\n" 
					<< "this is expected if generating MC files\n"
					<< "using default geometry from configuration file\n";

    return;
  }

  //......................................................................
  void Geometry::LoadGeometryFile(std::string gdmlfile, std::string rootfile)
  {
  
    struct stat sb;
    if (gdmlfile.empty() || stat(gdmlfile.c_str(), &sb)!=0)
      // Failed to resolve the file name
      throw cet::exception("NoGeometry") << "geometry file " << gdmlfile << " not found!";

    if (rootfile.empty() || stat(rootfile.c_str(), &sb)!=0)
      // Failed to resolve the file name
      throw cet::exception("NoGeometry") << "geometry file " << rootfile << " not found!";
 


    // clear the TPC array
    for (unsigned int i=0; i<fTPCs.size(); ++i) {
      if (fTPCs[i]) { delete fTPCs[i]; fTPCs[i] = 0; }
    }
    fTPCs.clear();

    // Open the GDML file, and convert it into ROOT TGeoManager
    // format.
    TGeoManager::Import(rootfile.c_str());

    std::vector<const TGeoNode*> path(8);
    path[0] = gGeoManager->GetTopNode();
    this->FindTPC(path, 0);

    std::sort(fTPCs.begin(), fTPCs.end(), tpc_sort);

    unsigned int chan = 0;
    for(unsigned int t = 0; t < fTPCs.size(); ++t){
      for(unsigned int tpcchan = 0; tpcchan < fTPCs[t]->Nchannels(); ++tpcchan){
	  TPCPlaneWireMap tpwm(t, tpcchan);
	  fChannelMap[chan] = tpwm;
	  ++chan;
      }// end loop over tpc channels
    }//end loops to make the channel map
  

    mf::LogWarning("LoadNewGeometry") << "New detector geometry loaded from "   
				      << "\n\t" << fROOTfile 
				      << "\n\t" << fGDMLfile << std::endl;

  }

  //......................................................................
  TGeoManager* Geometry::ROOTGeoManager() const
  {
    return gGeoManager;
  }
  
  //......................................................................
  unsigned int Geometry::Nchannels() const
  {
    return fChannelMap.size();
  }

  //......................................................................
  unsigned int Geometry::Nplanes(unsigned int tpc) const
  {
    if(tpc >= fTPCs.size() ){
      mf::LogWarning("TPCOutOfRange") << "requesting TPC " << tpc << " out of range, bail";
      assert(0);
    }

    return fTPCs[tpc]->Nplanes();
  }

  //......................................................................
  const std::string Geometry::GetWorldVolumeName() const
  {
    // For now, and possibly forever, this is a constant (given the
    // definition of "nodeNames" above).
    return std::string("volWorld");
  }

  //......................................................................
  const std::string Geometry::GetLArTPCVolumeName(unsigned int tpc) const
  {

    if(tpc >= fTPCs.size() ){
      mf::LogWarning("TPCOutOfRange") << "requesting TPC " << tpc << " out of range, bail";
      assert(0);
    }

    return std::string(fTPCs[tpc]->ActiveVolume()->GetName()); 
  }

  //......................................................................
  const std::string Geometry::GetCryostatVolumeName() const
  {
    // For now, and possibly forever, this is a constant (given the
    // definition of "nodeNames" above).
    return std::string("volCryostat");
  }

  //......................................................................
  double Geometry::DetHalfWidth(unsigned int tpc)  const 
  {
    return this->TPC(tpc).ActiveHalfWidth();
  }

  //......................................................................
  double Geometry::DetHalfHeight(unsigned int tpc) const 
  {
    return this->TPC(tpc).ActiveHalfHeight();
  }

  //......................................................................
  double Geometry::DetLength(unsigned int tpc) const
  { 
    return this->TPC(tpc).ActiveLength();
  }

  //......................................................................
  double Geometry::CryostatHalfWidth() const
  {
    // get the cryostat volume
    TGeoVolume *cv = gGeoManager->FindVolumeFast(this->GetCryostatVolumeName().c_str());
    return ((TGeoBBox*)cv->GetShape())->GetDX();
  }

  //......................................................................
  double Geometry::CryostatHalfHeight() const
  {
    // get the cryostat volume
    TGeoVolume *cv = gGeoManager->FindVolumeFast(this->GetCryostatVolumeName().c_str());
    return ((TGeoBBox*)cv->GetShape())->GetDY();
  }

  //......................................................................
  double Geometry::CryostatLength() const
  {
    // get the cryostat volume
    TGeoVolume *cv = gGeoManager->FindVolumeFast(this->GetCryostatVolumeName().c_str());
    return 2.0*((TGeoBBox*)cv->GetShape())->GetDZ();
  }

  //......................................................................
  // Boundaries of the cryostat in 3 pairs
  // [0]: -x
  // [1]: +x
  // [2]: -y
  // [3]: +y
  // [4]: -z
  // [5]: +z
  void Geometry::CryostatBoundaries(double* boundaries) const
  {
    TGeoVolume *v = gGeoManager->FindVolumeFast("volDetEnclosure");
    
    // loop over the nodes until we find the cryostat node
    for(int i = 0; i < v->GetNdaughters(); ++i){
      TGeoNode *n = v->GetNode(i);
      
      if(strncmp(n->GetName(), "volCryostat", strlen("volCryostat")) == 0){
	// get the half width, height, etc of the cryostat
	double halflength = ((TGeoBBox*)n->GetVolume()->GetShape())->GetDZ();
	double halfwidth  = ((TGeoBBox*)n->GetVolume()->GetShape())->GetDX();
	double halfheight = ((TGeoBBox*)n->GetVolume()->GetShape())->GetDY();

	double posW[3] = {0.};
	double negW[3] = {0.};
	double pos[3] = {halfwidth, halfheight, halflength};
	double neg[3] = {-halfwidth, -halfheight, -halflength};
        
	n->LocalToMaster(pos, posW);
	n->LocalToMaster(neg, negW);
	
	boundaries[0] = negW[0];
	boundaries[1] = posW[0];
	boundaries[2] = negW[1];
	boundaries[3] = posW[1];
	boundaries[4] = negW[2];
	boundaries[5] = posW[2];

	break;
      }
    }// end loop over daughter nodes of volDetEnclosure

    return;
  }

  //......................................................................
  // This method returns the distance between the specified planes.
  // p1 < p2
  double Geometry::PlanePitch(unsigned int p1, 
			      unsigned int p2, 
			      unsigned int tpc) const
  { 
    return this->TPC(tpc).PlanePitch(p1, p2);
  }
      
  //......................................................................
  // This method returns the distance between the specified wires.
  // w1 < w2.  The wires are assumed to be on the same plane
  double Geometry::WirePitch(unsigned int w1,  
			     unsigned int w2,  
			     unsigned int plane,
			     unsigned int tpc) const
  { 
    return this->TPC(tpc).WirePitch(w1,w2,plane);    
  }

  //......................................................................
  //
  // Return the ranges of x,y and z for the "world volume" that the
  // entire geometry lives in. If any pointers are 0, then those
  // coordinates are ignored.
  //
  // \param xlo : On return, lower bound on x positions
  // \param xhi : On return, upper bound on x positions
  // \param ylo : On return, lower bound on y positions
  // \param yhi : On return, upper bound on y positions
  // \param zlo : On return, lower bound on z positions
  // \param zhi : On return, upper bound on z positions
  //
  void Geometry::WorldBox(double* xlo, double* xhi,
			  double* ylo, double* yhi,
			  double* zlo, double* zhi) const
  {
    const TGeoShape* s = gGeoManager->GetVolume("volWorld")->GetShape();
    assert(s);

    if (xlo || xhi) {
      double x1, x2;
      s->GetAxisRange(1,x1,x2); if (xlo) *xlo = x1; if (xhi) *xhi = x2;
    }
    if (ylo || yhi) {
      double y1, y2;
      s->GetAxisRange(2,y1,y2); if (ylo) *ylo = y1; if (yhi) *yhi = y2;
    }
    if (zlo || zhi) {
      double z1, z2;
      s->GetAxisRange(3,z1,z2); if (zlo) *zlo = z1; if (zhi) *zhi = z2;
    }
  }

  //......................................................................
  const TVector3 Geometry::GetTPCFrontFaceCenter() const
  {
    return TVector3( 0.5 * this->DetHalfWidth(), 0 , 0 );
  }

  //......................................................................
  const std::string Geometry::VolumeName(TVector3 point)
  {
    // check that the given point is in the World volume at least
    TGeoVolume *volWorld = gGeoManager->FindVolumeFast(this->GetWorldVolumeName().c_str());
    double halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
    double halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
    double halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
    if(TMath::Abs(point.x()) > halfwidth  ||
       TMath::Abs(point.y()) > halfheight ||
       TMath::Abs(point.z()) > halflength
       ){
      if(fNBadVol%(int)(1e6)==0){
	mf::LogWarning("GeometryBadInputPoint") << "point (" << point.x() << ","
						<< point.y() << "," << point.z() << ") "
						<< "is not inside the world volume "
						<< " half width = " << halfwidth
						<< " half height = " << halfheight
						<< " half length = " << halflength
						<< " returning unknown volume name";
      }
      fNBadVol++;
      const std::string unknown("unknownVolume");
      return unknown;
    }

    const std::string name(gGeoManager->FindNode(point.x(), point.y(), point.z())->GetVolume()->GetName());
    return name;
  }

  //......................................................................
  const std::string Geometry::MaterialName(TVector3 point)
  {
    // check that the given point is in the World volume at least
    TGeoVolume *volWorld = gGeoManager->FindVolumeFast(this->GetWorldVolumeName().c_str());
    double halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
    double halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
    double halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
    if(TMath::Abs(point.x()) > halfwidth  ||
       TMath::Abs(point.y()) > halfheight ||
       TMath::Abs(point.z()) > halflength
       ){ 
      if (fNBadVol%(int)(1e6)==0){
	    mf::LogWarning("GeometryBadInputPoint") << "point (" << point.x() << ","
						    << point.y() << "," << point.z() << ") "
						    << "is not inside the world volume "
						    << " half width = " << halfwidth
						    << " half height = " << halfheight
						    << " half length = " << halflength
						    << " returning unknown material name";
      }
      const std::string unknown("unknownMaterial");
      return unknown;
    }
    
    const std::string name(gGeoManager->FindNode(point.x(), point.y(), point.z())->GetMedium()->GetMaterial()->GetName());
    return name;
  }

  //......................................................................
  void Geometry::FindTPC(std::vector<const TGeoNode*>& path,
			 unsigned int depth)
  {

    const char* nm = path[depth]->GetName();
    if( (strncmp(nm, "volTPC", 6) == 0) ){
      this->MakeTPC(path,depth);
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
      this->FindTPC(path, deeper);
    }
  
  }

  //......................................................................
  void Geometry::MakeTPC(std::vector<const TGeoNode*>& path, int depth) 
  {
    fTPCs.push_back(new TPCGeo(path, depth));
  }

  //......................................................................
  //
  // Return the geometry description of the ith plane in the detector.
  //
  // \param tpc : input plane number, starting from 0
  // \returns plane geometry for ith plane
  //
  // \throws geo::Exception if "tpc" is outside allowed range
  //
  const TPCGeo& Geometry::TPC(unsigned int tpc) const 
  {
    if(tpc < 0 || tpc >= fTPCs.size()) 
      throw cet::exception("TPCOutOfRange") << "tpc " << tpc << " does not exist";

    return *fTPCs[tpc];
  }

  //......................................................................
  const TPCGeo& Geometry::PositionToTPC(double *worldLoc) const
  {
    // boundaries of the TPC in the world volume are organized as
    // [0] = -x
    // [1] = +x
    // [2] = -y
    // [3] = +y
    // [4] = -z
    // [5] = +z
    static std::vector<double> tpcBoundaries(this->NTPC()*6);

    bool firstCalculation = true;

    if ( firstCalculation ){
      firstCalculation = false;
      double origin[3] = {0.};
      double world[3] = {0.};
      for(unsigned int t = 0; t < this->NTPC(); ++t){
	this->TPC(t).LocalToWorld(origin, world);
	// y and z values are easy and can be figured out using the TPC origin
	// the x values are a bit trickier, at least the -x value seems to be
	tpcBoundaries[0+t*6] =  world[0] - this->TPC(t).HalfWidth();
	tpcBoundaries[1+t*6] =  world[0] + this->TPC(t).HalfWidth();
	tpcBoundaries[2+t*6] =  world[1] - this->TPC(t).HalfHeight();
	tpcBoundaries[3+t*6] =  world[1] + this->TPC(t).HalfHeight();
	tpcBoundaries[4+t*6] =  world[2] - 0.5*this->TPC(t).Length();
	tpcBoundaries[5+t*6] =  world[2] + 0.5*this->TPC(t).Length();
      }
    }// end if this is the first calculation

    // locate the desired TPC
    unsigned int tpc = UINT_MAX;
    for(unsigned int t = 0; t < this->NTPC(); ++t){
      if(worldLoc[0] >= tpcBoundaries[0+t*6] &&
	 worldLoc[0] <= tpcBoundaries[1+t*6] && 
	 worldLoc[1] >= tpcBoundaries[2+t*6] && 
	 worldLoc[1] <= tpcBoundaries[3+t*6] && 
	 worldLoc[2] >= tpcBoundaries[4+t*6] && 
	 worldLoc[2] <= tpcBoundaries[5+t*6] ){
	tpc = t;
	break;
      }
    }

    if(tpc == UINT_MAX)
      throw cet::exception("Can't find TPC") << "for position (" << worldLoc[0] << ","
					     << worldLoc[1] << "," << worldLoc[2] << ")";
			
    return this->TPC(tpc);
  }

  //......................................................................
  //
  // Return the geometry description of the ith plane in the detector.
  //
  // \param p : input plane number, starting from 0
  // \returns plane geometry for ith plane
  //
  // \throws geo::Exception if "i" is outside allowed range
  //
  const PlaneGeo& Geometry::Plane(unsigned int p, unsigned int tpc) const 
  {
    if(tpc < 0 || tpc >= fTPCs.size()) 
      throw cet::exception("TPCOutOfRange") << "tpc " << tpc << " does not exist";

    return fTPCs[tpc]->Plane(p);
  }

  //......................................................................
  unsigned int Geometry::Nwires(unsigned int p, unsigned int tpc) const 
  {
    if(tpc >= fTPCs.size()) 
      throw cet::exception("PlaneOutOfRange") << "plane " << p << " does not exist";

    return fTPCs[tpc]->Plane(p).Nwires();
  }

  //......................................................................
  //
  // Return the total mass of the detector
  //
  //
  double Geometry::TotalMass(const char *vol) const
  {
    //the TGeoNode::GetVolume() returns the TGeoVolume of the detector outline
    //and ROOT calculates the mass in kg for you
    TGeoVolume *gvol = gGeoManager->FindVolumeFast(vol);
    if(gvol) return gvol->Weight();

    throw cet::exception("NoVolumeFound") << "could not find specified volume " << vol 
					  << " returning weight of active volume" << std::endl;
    

  }

  //......................................................................
  //
  // Return the column density between 2 points
  //
  // \param p1  : pointer to array holding xyz of first point in world coordinates
  // \param p2  : pointer to array holding xyz of second point in world coorinates
  //
  double Geometry::MassBetweenPoints(double *p1, double *p2) const
  {

    //The purpose of this method is to determine the column density
    //between the two points given.  Do that by starting at p1 and 
    //stepping until you get to the node of p2.  calculate the distance
    //between the point just inside that node and p2 to get the last
    //bit of column density
    double columnD = 0.;

    //first initialize a track - get the direction cosines
    double length = TMath::Sqrt(TMath::Power(p2[0]-p1[0], 2.)
				+ TMath::Power(p2[1]-p1[1], 2.)
				+ TMath::Power(p2[2]-p1[2], 2.));
    double dxyz[3] = {(p2[0]-p1[0])/length, (p2[1]-p1[1])/length, (p2[2]-p1[2])/length}; 

    gGeoManager->InitTrack(p1,dxyz);

    //might be helpful to have a point to a TGeoNode
    TGeoNode *node = gGeoManager->GetCurrentNode();

    //check that the points are not in the same volume already.  
    //if they are in different volumes, keep stepping until you 
    //are in the same volume as the second point
    while(!gGeoManager->IsSameLocation(p2[0], p2[1], p2[2])){
      gGeoManager->FindNextBoundary();
      columnD += gGeoManager->GetStep()*node->GetMedium()->GetMaterial()->GetDensity();
    
      //the act of stepping puts you in the next node and returns that node
      node = gGeoManager->Step();
    }//end loop to get to volume of second point

    //now you are in the same volume as the last point, but not at that point.
    //get the distance between the current point and the last one
    const double *current = gGeoManager->GetCurrentPoint();
    length = TMath::Sqrt(TMath::Power(p2[0]-current[0], 2.)
			 + TMath::Power(p2[1]-current[1], 2.)
			 + TMath::Power(p2[2]-current[2], 2.));
    columnD += length*node->GetMedium()->GetMaterial()->GetDensity();

    return columnD;
  }

  //......................................................................
  const WireGeo& Geometry::ChannelToWire(unsigned int  channel,
					 unsigned int &tpc, 
					 unsigned int &plane,
					 unsigned int &wire) const
  {
    //the conversion from channel number to plane and wire is straightforward
    ChannelMap::const_iterator itr = fChannelMap.find(channel);
    if( itr == fChannelMap.end() ){
      throw cet::exception("geo::ChannelToWire(): BAD CHANNEL Number ") << "channel " << channel << " "
									<< " not found in map";
      
    }

    tpc   = itr->second.first;
    unsigned int tpcchan = itr->second.second;

    return this->TPC(tpc).ChannelToWire(tpcchan, plane, wire);
  }


  //......................................................................
  void Geometry::ChannelToWireFast(unsigned int  channel,
				   unsigned int &tpc, 
				   unsigned int &plane,
				   unsigned int &wire) const
  {
    //Analagous with PlaneWireToChannelFast - lookup assuming a heirachical
    // wire labelling scheme (see comments above that method) - Ben J Oct 2011

    static std::vector<std::vector<unsigned int> > FirstChannelInNextPlane;  // we gain efficiency by storing both of these
    static std::vector<std::vector<unsigned int> > FirstChannelInThisPlane;  //  their contents are self explanatory
   
    static unsigned int TopChannel = 0;
    static unsigned int ntpc = this->NTPC();
    
    // run intialization on first call only
    static bool FirstCall = true;
    if(FirstCall){
      
      FirstCall=false;
      
      FirstChannelInThisPlane.resize(ntpc);
      FirstChannelInNextPlane.resize(ntpc);
      
      for(unsigned int TPCCount = 0; TPCCount != ntpc; ++TPCCount)
	for(unsigned int PCount = 0; PCount != TPC(TPCCount).Nplanes(); ++PCount){
	  unsigned int WiresThisPlane = TPC(TPCCount).Plane(PCount).Nwires();
	  FirstChannelInThisPlane.at(TPCCount).push_back(TopChannel);
	  TopChannel += WiresThisPlane;
	  FirstChannelInNextPlane.at(TPCCount).push_back(TopChannel);
	}

    } // end of initialization


    // This part runs every time

    // first check if this channel ID is legal
    if((channel<0)||(channel>TopChannel)){
      throw cet::exception("Geometry::ILLEGAL CHANNEL ID ")<<"for channel" << channel	<< "\n returning UINT_MAX";
      wire = UINT_MAX;
      plane = UINT_MAX;
      tpc = UINT_MAX;
    }
    
    // then go find which plane and tpc it is in from the information we stored earlier
    for(unsigned int tpcloop = 0; tpcloop != ntpc; ++tpcloop)
      for(unsigned int planeloop = 0; planeloop != FirstChannelInNextPlane.at(tpcloop).size(); ++planeloop){
	if(channel<FirstChannelInNextPlane[tpcloop][planeloop]){
	  tpc   = tpcloop;
	  plane = planeloop;
	  wire  = channel - FirstChannelInThisPlane[tpcloop][planeloop];
	  return;
	}	    
      }
    
  }

  //......................................................................
  unsigned int Geometry::PlaneWireToChannel(unsigned int plane,
					    unsigned int wire,
					    unsigned int tpc) const
  {
    // first get the channel from the TPC
    unsigned int tpcchan = this->TPC(tpc).PlaneWireToTPCChannel(plane, wire);

    TPCPlaneWireMap tpwm(tpc,tpcchan); 

    ChannelMap::const_iterator itr = fChannelMap.begin();
    while(itr != fChannelMap.end() ){
      if(itr->second == tpwm) return itr->first;
      itr++;
    }

    // if we got here the tpc, plane, and wire set are not in the geometry

    throw cet::exception("Geometry::PlaneWireToChannel: NO CHANNEL FOUND ") << "for tpc,plane,wire: "
									    << tpc << "," << plane
									    << "," << wire
									    << "\n returning UINT_MAX";

    return UINT_MAX;

  }



  //---------------------------------------

  // The two routines below are attempts to speed up the simulation
  //  by avoiding computationally intensive geometry calculations.
  //  The results are valid assuming the wireplanes have comprised of
  //  straight, parallel wires with constant pitch across the entire plane,
  //  with a heirachical numbering scheme -  Ben J Oct 2011
  unsigned int Geometry::NearestChannelFast(double *worldPos, 
					    unsigned int PlaneNo, 
					    unsigned int TPCNo) const
  {

    // All these variables are static - avoid looking up too many things on each
    //  call to the routine

    static bool FirstCalc=true;    
    static unsigned int ntpc;  // number of TPCs - store statically so we dont have to keep looking it up

    // Distance (0,0,0) to first wire along orth vector per plane per TPC
    //    (see further notes below)
    static std::vector< std::vector<float> >   FirstWireProj;

    // Unit vectors orthogonal to wires in each plane - stored as 2 components
    //    to avoid having to invoke any bulky TObjects / CLHEP vectors etc
    //    (see further notes below)
    static std::vector< std::vector<float > > OrthVectorsY;
    static std::vector< std::vector<float > > OrthVectorsZ;

    // Number of wires in each plane - for range checking after calculation
    static std::vector< std::vector<float> > WireCounts;

    // Initialization - run only on the first call
    if(FirstCalc){

      mf::LogInfo("Geometry") <<"NearestChannelFast initializing...";
      
      FirstCalc=false;
      
      ntpc = NTPC();
      
      // Size up all the vectors 
      WireCounts.resize(ntpc);
      FirstWireProj.resize(ntpc);
      OrthVectorsY.resize(ntpc);
      OrthVectorsZ.resize(ntpc);
      
      for(unsigned int TPCCount = 0; TPCCount != ntpc; ++TPCCount){
	unsigned int PlanesThisTPC = Nplanes(TPCCount);
	WireCounts   .at( TPCCount ).resize(PlanesThisTPC);
	FirstWireProj.at( TPCCount ).resize(PlanesThisTPC);
	OrthVectorsY .at( TPCCount ).resize(PlanesThisTPC);
	OrthVectorsZ .at( TPCCount ).resize(PlanesThisTPC);
      }
	
      for(unsigned int TPCCount = 0; TPCCount != NTPC(); ++TPCCount)
	for(unsigned int PlaneCount = 0; PlaneCount != Nplanes(TPCCount); ++PlaneCount){
	  double ThisWirePitch = WirePitch(0,1,PlaneCount,TPCCount);
	  WireCounts.at(TPCCount).at(PlaneCount) = Nwires(PlaneCount, TPCCount);
	  
	  double  WireCentre1[3] = {0.,0.,0.};
	  double  WireCentre2[3] = {0.,0.,0.};
	  
	  double th = this->TPC(TPCCount).Plane(PlaneCount).Wire(0).ThetaZ();
	  double sth = sin(th), cth = cos(th);
	  
	  this->TPC(TPCCount).Plane(PlaneCount).Wire(0).GetCenter(WireCentre1,0);
	  this->TPC(TPCCount).Plane(PlaneCount).Wire(1).GetCenter(WireCentre2,0);
	  
	  // figure out if we need to flip the orthogonal vector 
	  // (should point from wire n -> n+1)
	  double OrthY=cth, OrthZ=-sth;
	  if(((WireCentre2[1]-WireCentre1[1])*OrthY 
	      + (WireCentre2[2]-WireCentre1[2])*OrthZ) < 0){
	    OrthZ *= -1;
	    OrthY *= -1;
	  }
	  
	  // Overall we are trying to build an expression that looks like
	  //       int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);     
	  // That runs as fast as humanly possible.
	  //
	  // Casting to an int is much faster than c rounding commands like floor().  We have to add 0.5
	  // to account for rounding up as well as down.  floor(A) ~ (int)(A+0.5).  We account for the
	  // 0.5 in the first wire constant to avoid adding it every time.
	  //
	  // We can also predivide everything by the wire pitch so we don't do this in the loop
	  //
	  // Putting this together into the useful constants we will use later per plane and tpc:
	  OrthVectorsY.at( TPCCount).at( PlaneCount ) = OrthY / ThisWirePitch;
	  OrthVectorsZ.at( TPCCount).at( PlaneCount ) = OrthZ / ThisWirePitch;
	  
	  FirstWireProj.at( TPCCount).at( PlaneCount) = (WireCentre1[1] * OrthY + WireCentre1[2] * OrthZ) / ThisWirePitch - 0.5;
	}
    }

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
     
    int NearestWireNumber = (int)(worldPos[1]*OrthVectorsY[TPCNo][PlaneNo] 
				  + worldPos[2]*OrthVectorsZ[TPCNo][PlaneNo]      
				  - FirstWireProj[TPCNo][PlaneNo]);

    // If we are outside of the wireplane range, throw an exception
    //  (this response maintains consistency with the previous
    //   implementation based on geometry lookup)
    
    if((NearestWireNumber < 0)||(NearestWireNumber >= WireCounts[TPCNo][PlaneNo])){
      throw cet::exception("Can't Find Nearest Channel") << "for position (" 
							 << worldPos[0] << ","
							 << worldPos[1] << "," 
							 << worldPos[2] << ")"
							 << "\n return UINT_MAX";
      return UINT_MAX;
    }
    
    // This method is supposed to return a channel number rather than
    //  a wire number.  Perform the conversion here (although, maybe
    //  faster if we deal in wire numbers rather than channel numbers?)

    return PlaneWireToChannelFast(PlaneNo, NearestWireNumber, TPCNo);
  }

  //--------------------------------------
  // This method returns the channel number, assuming the numbering scheme
  // is heirachical - that is, channel numbers run in order, for example:
  //                                             (Ben J Oct 2011)                   
  //                    Wire1     | 0
  //           Plane1 { Wire2     | 1
  //    TPC1 {          Wire3     | 2
  //           Plane2 { Wire1     | 3   increasing channel number
  //                    Wire2     | 4     (with no gaps)
  //    TPC2 { Plane1 { Wire1     | 5
  //           Plane2 { Wire1     | 6
  //                    Wire2     v 7
  //
 
  unsigned int Geometry::PlaneWireToChannelFast(unsigned int plane,
						unsigned int wire,
						unsigned int tpc) const
  {
    // Store these statically so we don't have to keep looking up / calculating

    // The number of wires in all the tpcs and planes up to this one in the heirachy
    static std::vector<std::vector<unsigned int> > PlaneBaselines;

    // The number of wires in this plane in the heirachy
    static std::vector<std::vector<unsigned int> > WiresPerPlane;
    
    static unsigned int ntpc=NTPC();
    
    // run intialization on first call only
    static bool FirstCall = true;
    if(FirstCall){
      
      FirstCall=false;
      
      PlaneBaselines.resize(ntpc);
      WiresPerPlane.resize(ntpc);
      
      int RunningTotal=0;
      
      for(unsigned int TPCCount=0; TPCCount!=ntpc; ++TPCCount)
	for(unsigned int PCount=0; PCount!=TPC(TPCCount).Nplanes(); ++PCount){
	  int WiresThisPlane = TPC(TPCCount).Plane(PCount).Nwires();	
	  WiresPerPlane.at(TPCCount).push_back(WiresThisPlane);
	  PlaneBaselines.at(TPCCount).push_back(RunningTotal);
	  
	  RunningTotal+=WiresThisPlane;
	}
    }

    // This is the actual lookup part - first make sure coordinates are legal
    if((tpc < ntpc) &&
       (plane < WiresPerPlane[tpc].size()) &&
       (wire < WiresPerPlane[tpc][plane])){
      // if the channel has legal coordinates, its ID is given by the wire
      //   number above the number of wires in lower planes and lowers tpcs
      return PlaneBaselines[tpc][plane] + wire;
    }
    else{  // if the coordinates were bad, throw an exception
     
      throw cet::exception("Geometry::PlaneWireToChannelFast: NO CHANNEL FOUND ") << "for tpc,plane,wire: "
										  << tpc << "," << plane
										  << "," << wire
										  << "\n returning UINT_MAX";
      return UINT_MAX;
    }
  }

  //......................................................................
  unsigned int Geometry::NearestChannel(double* worldPos) const
  {
    // This routine will throw a cet::Exception if it cannot locate
    // the nearest channel, best to put try/catch around calls to it

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
    static unsigned int nearest = UINT_MAX;
    static double closeEnough = 0.;
    static bool firstCalculation = true;


    if ( firstCalculation ){
      firstCalculation = false;
      // What do we mean by "close enough" to the previous point?
      // For now, let's take it to be one-tenth the of the wire
      // spacing; this is typically the same size as the voxels used
      // in the simulation.
      closeEnough = this->WirePitch( 0, 1, 0, 0 ) / 10.;
    }
    else{
      if ( fabs(worldPos[0] - lastPos0) < closeEnough  &&
	   fabs(worldPos[1] - lastPos1) < closeEnough  &&
	   fabs(worldPos[2] - lastPos2) < closeEnough  ){ 
	return nearest; 
      }
    }

    // If we get to this line, we're doing the full calculation. Save
    // the current point.
    lastPos0 = worldPos[0];
    lastPos1 = worldPos[1];
    lastPos2 = worldPos[2];

    try{
      const TPCGeo &tpc = this->PositionToTPC(worldPos);
      nearest = tpc.NearestChannel(worldPos);
      return nearest;
    }
    catch(cet::exception &e){
      mf::LogWarning("FindTPCException") << e;

      // if we got here then we didn't find the closest channel so
      // return something else
      throw cet::exception("Can't Find Nearest Channel") << "for position (" 
							 << worldPos[0] << ","
							 << worldPos[1] << "," 
							 << worldPos[2] << ")"
							 << "\n return UINT_MAX";
    }

    return UINT_MAX;
  }

  // Functions to allow determination if two wires intersect, and if so where.
  // This is useful information during 3D reconstruction.
  //......................................................................
  bool Geometry::ValueInRange(double value, double min, double max)
  {
     if(min>max) std::swap(min,max);//protect against funny business due to wire angles
     return (value>=min) && (value<=max);
  }

  //......................................................................
  void Geometry::WireEndPoints(unsigned int tpc,
			       unsigned int plane, 
			       unsigned int wire, 
			       double *xyzStart, 
			       double *xyzEnd)
  {  
    double halfL = this->TPC(tpc).Plane(plane).Wire(wire).HalfL();//half-length of wire
    this->TPC(tpc).Plane(plane).Wire(wire).GetCenter(xyzStart,halfL);
    this->TPC(tpc).Plane(plane).Wire(wire).GetCenter(xyzEnd,-1.0*halfL);
    
    if(xyzEnd[2]<xyzStart[2]){
      //ensure that "End" has higher z-value than "Start"
      std::swap(xyzStart[0],xyzEnd[0]);
      std::swap(xyzStart[1],xyzEnd[1]);
      std::swap(xyzStart[2],xyzEnd[2]);
    }
    if(xyzEnd[1]<xyzStart[1] && fabs(xyzEnd[2]-xyzStart[2])<0.01){
      // if wire is vertical ensure that "End" has higher y-value than "Start"
      std::swap(xyzStart[0],xyzEnd[0]);
      std::swap(xyzStart[1],xyzEnd[1]);
      std::swap(xyzStart[2],xyzEnd[2]);
    }
    
     return;
  }
   
  //......................................................................
  bool Geometry::ChannelsIntersect(unsigned short c1, 
				   unsigned short c2, 
				   double &y, 
				   double &z)
  {
    unsigned int tpc1,plane1,wire1;
    unsigned int tpc2,plane2,wire2;
    this->ChannelToWire(c1,tpc1,plane1,wire1);
    this->ChannelToWire(c2,tpc2,plane2,wire2);

    if(tpc1 != tpc2){
      mf::LogWarning("ChannelsIntersect") << "attempting to find intersection between wires"
					  << " from different TPCs, return false";
      return false;
    }

    double wire1_Start[3] = {0.};
    double wire1_End[3]   = {0.};
    double wire2_Start[3] = {0.};
    double wire2_End[3]   = {0.};
	
    this->WireEndPoints(tpc1,plane1,wire1,wire1_Start,wire1_End);
    this->WireEndPoints(tpc1,plane2,wire2,wire2_Start,wire2_End);

    if(plane1==plane2){
      mf::LogWarning("ChannelsIntersect") << "You are comparing two wires in the same plane!";
      return false;
    }

    // if endpoint of one input wire is within range of other input wire in 
    // BOTH y AND z, wires overlap 
    bool overlapY = this->ValueInRange(wire1_Start[1],wire2_Start[1],wire2_End[1]) ||
      this->ValueInRange(wire1_End[1],wire2_Start[1],wire2_End[1]);
    
    bool overlapZ = this->ValueInRange(wire1_Start[2],wire2_Start[2],wire2_End[2]) ||
      this->ValueInRange(wire1_End[2],wire2_Start[2],wire2_End[2]);
    
    // reverse ordering of wires...this is necessitated for now due to precision 
    // of placement of wires 
    bool overlapY_reverse = this->ValueInRange(wire2_Start[1],wire1_Start[1],wire1_End[1]) ||
      this->ValueInRange(wire2_End[1],wire1_Start[1],wire1_End[1]);
    
    bool overlapZ_reverse = this->ValueInRange(wire2_Start[2],wire1_Start[2],wire1_End[2]) ||
      this->ValueInRange(wire2_End[2],wire1_Start[2],wire1_End[2]);
 
    // override y overlap checks if a vertical plane exists:
    if( (this->TPC(tpc1).Plane(plane1).Wire(wire1).ThetaZ() == M_PI/2) 
	|| (this->TPC(tpc1).Plane(plane2).Wire(wire2).ThetaZ() == M_PI/2)){
      overlapY=true;	
      overlapY_reverse=true;
    }

    //catch to get vertical wires, where the standard overlap might not work, Andrzej
    if(fabs(wire2_Start[2]-wire2_End[2])<0.01) overlapZ=overlapZ_reverse;

    if(overlapY && overlapZ){
      this->IntersectionPoint(wire1, wire2, 
			      plane1, plane2, 
			      tpc1,
			      wire1_Start, wire1_End, 
			      wire2_Start, wire2_End, 
			      y, z);
      return true;
    }
    else if(overlapY_reverse && overlapZ_reverse){
      this->IntersectionPoint(wire2, wire1, 
			      plane2, plane1, 
			      tpc1,
			      wire2_Start, wire2_End, 
			      wire1_Start, wire1_End, 
			      y, z);
      return true;
    }
    else return false;
    
  }
   
  //......................................................................
  // This function is called if it is determined that two wires in a single TPC must overlap.
  // To determine the yz coordinate of the wire intersection, we need to know the 
  // endpoints of both wires in xyz-space, and also their orientation (angle), and the 
  // inner dimensions of the TPC frame.
  // Note: This calculation is entirely dependent  on an accurate GDML description of the TPC!
  // Mitch - Feb., 2011

  void Geometry::IntersectionPoint(unsigned int wire1, 
				   unsigned int wire2, 
				   unsigned int plane1, 
				   unsigned int plane2,
				   unsigned int tpc,
                                   double start_w1[3], 
				   double end_w1[3], 
				   double start_w2[3], 
				   double end_w2[3], 
                                   double &y, double &z)
  {

    //angle of wire1 wrt z-axis in Y-Z plane...in radians
    double angle1 = this->TPC(tpc).Plane(plane1).Wire(wire1).ThetaZ();
    //angle of wire2 wrt z-axis in Y-Z plane...in radians
    double angle2 = this->TPC(tpc).Plane(plane2).Wire(wire2).ThetaZ();
    
    if(angle1==angle2) return;//comparing two wires in the same plane...pointless.

    //coordinates of "upper" endpoints...(z1,y1) = (a,b) and (z2,y2) = (c,d) 
    double a,b,c,d;
    double angle,anglex;
    
    // below is a special case of calculation when one of the planes is vertical. 
    angle1<angle2 ? angle = angle1 : angle = angle2;//get angle closest to the z-axis
    
    // special case, one plane is vertical
    if(angle1 == M_PI/2 || angle2 == M_PI/2){
      mf::LogInfo("Geometry") << "calculating for vertical plane ";
      if(angle1 == M_PI/2){
		
	anglex=(angle2-M_PI/2);
	a=end_w1[2];
	b=end_w1[1];
	c=end_w2[2];
	d=end_w2[1];
	// the if below can in principle be replaced by the sign of anglex (inverted) 
	// in the formula for y below. But until the geometry is fully symmetric in y I'm 
	// leaving it like this. Andrzej
	if((anglex) > 0 ) b=start_w1[1];
		    
      }
      else if(angle2==M_PI/2){
	anglex=(angle1-M_PI/2);
	a=end_w2[2];
	b=end_w2[1];
	c=end_w1[2];
	d=end_w1[1];
	// the if below can in principle be replaced by the sign of anglex (inverted) 
	// in the formula for y below. But until the geometry is fully symmetric in y I'm 
	// leaving it like this. Andrzej
	if((anglex) > 0 ) b=start_w2[1];  
      }

      y=b+((c-a) - (b-d)*tan(anglex))/tan(anglex);
      z=a;   // z is defined by the wire in the vertical plane
      
      return;
    }

    // end of vertical case
   
    z=0;y=0;
                                                                      
    if(angle1<(TMath::Pi()/2.0)){
      c=end_w1[2];
      d=end_w1[1];
      a=start_w2[2];
      b=start_w2[1];
    }
    else{
      c=end_w2[2];
      d=end_w2[1];
      a=start_w1[2];
      b=start_w1[1];
    }
    
    //Intersection point of two wires in the yz plane is completely
    //determined by wire endpoints and angle of inclination.
    z = 0.5 * ( c + a + (b-d)/TMath::Tan(angle) );
    y = 0.5 * ( b + d + (a-c)*TMath::Tan(angle) );
    
    return;

  }
    
  // Added shorthand function where start and endpoints are looked up automatically
  //  - whether to use this or the full function depends on optimization of your
  //    particular algorithm.  Ben J, Oct 2011
  //--------------------------------------------------------------------
  void Geometry::IntersectionPoint(unsigned int wire1, 
				   unsigned int wire2, 
				   unsigned int plane1, 
				   unsigned int plane2,
				   unsigned int tpc, 
                                   double &y, double &z)
  {
    double WireStart1[3], WireStart2[3], WireEnd1[3], WireEnd2[3];
    WireEndPoints(tpc, plane1, wire1, WireStart1, WireEnd1);
    WireEndPoints(tpc, plane2, wire2, WireStart2, WireEnd2);
    IntersectionPoint(wire1, wire2, plane1, plane2, tpc,
		      WireStart1, WireEnd1, WireStart2, WireEnd2, y, z);		     
  }

  //......................................................................
  bool Geometry::CrossesVol( double xyz[],
			     double dxyz[],
			     double point[] )
  {
    
    double extra = 5000;
    double height = 0;
    
    double x_lo = 0                       - extra;
    double x_hi = 2*this->DetHalfWidth()  + extra;
    double y_lo = - this->DetHalfHeight() - extra + height;
    double y_hi =   this->DetHalfHeight() + extra + height;
    double z_lo = 0                       - extra;
    double z_hi =   this->DetLength()     + extra;
    
    //std::cout << "  Box size (x0,x1;y0,y1;z0,z1) = (" << x_lo << "," << x_hi << ";" << y_lo << "," << y_hi << ";" << z_lo << "," << z_hi << ")" << std::endl;
    return CrossesBoundary(xyz, dxyz, x_lo, x_hi, y_lo, y_hi, z_lo, z_hi, point); 
  }
  
} // namespace geo
