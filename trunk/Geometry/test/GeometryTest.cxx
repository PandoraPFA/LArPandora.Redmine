////////////////////////////////////////////////////////////////////////
/// \file    GeometryTest.cxx
/// \brief   unit tests for the Geometry service
/// \version $Id: GeometryTest.cxx,v 1.1 2011/02/17 01:45:48 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <vector>

// ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TVector3.h"
#include "TNtuple.h"
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TMath.h"

// LArSoft includes
#include "Geometry/test/GeometryTest.h"
#include "Geometry/geo.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


namespace geo{

  //......................................................................
  GeometryTest::GeometryTest(fhicl::ParameterSet const& pset) 
  {
  }

  //......................................................................
  GeometryTest::~GeometryTest()
  {
  }

  //......................................................................
  void GeometryTest::analyze(art::Event const& evt)
  {
    art::ServiceHandle<geo::Geometry> geom;
    
    try{
      std::cout << "Wire Rmax  "  << geom->Plane(1).Wire(10).RMax()     << std::endl;
      std::cout << "Wire length " << 2.*geom->Plane(1).Wire(10).HalfL() << std::endl;
      std::cout << "Wire Rmin  "  << geom->Plane(1).Wire(10).RMin()     << std::endl;
      std::cout << "Total mass "  << geom->TotalMass()                  << std::endl;
      std::cout << "Cryostat dimensions: " << 2.*geom->CryostatHalfWidth() 
		<< " x "          << 2.*geom->CryostatHalfHeight()
		<< " x "          << geom->CryostatLength()             << std::endl;

      double cryobound[6] = {0.};
      geom->CryostatBoundaries(cryobound);
      std::cout << "Cryostat boundaries are at:\n"
		<< "\t-x:" << cryobound[0] << " +x:" << cryobound[1]
		<< "\t-y:" << cryobound[2] << " +y:" << cryobound[3]
		<< "\t-z:" << cryobound[4] << " +z:" << cryobound[5] << std::endl;

      std::cout << "test TPC methods ..." << std::endl;
      testTPC();
      std::cout << "complete." << std::endl;

      std::cout << "test channel to plane wire and back ..." << std::endl;
      testChannelToWire();
      std::cout << "complete." << std::endl;

      std::cout << "test find plane centers..." << std::endl;
      testFindPlaneCenters();
      std::cout << "complete." << std::endl;

      std::cout << "testProject...";
      testProject();
      std::cout << "complete." << std::endl;

      std::cout << "testWirePos...";
      testWirePos();
      std::cout << "complete." << std::endl;
      
      std::cout << "testWirePitch...";
      testWirePitch();
      std::cout << "complete." << std::endl;

      std::cout << "testPlanePitch...";
      testPlanePitch();
      std::cout << "complete." << std::endl;

      std::cout << "testStepping...";
      testStepping();
      std::cout << "complete." << std::endl;

    }
    catch (...) {
      abort();
    }
    
    return;
  }

  //......................................................................
  void GeometryTest::testTPC()
  {
    art::ServiceHandle<geo::Geometry> geom;

    std::cout << "\tThere are " << geom->NTPC() << " TPCs in the detector" << std::endl;
    
    for(size_t t = 0; t < geom->NTPC(); ++t){
      std::cout << std::endl << "\t\tTPC " << t << " has " 
		<< geom->Nplanes(t) << " planes." << std::endl;
      for(size_t p = 0; p < geom->Nplanes(t); ++p)
	std::cout << std::endl << "\t\tPlane " << p << " has " 
		  << geom->Plane(p).Nwires() << " wires and is at (x,y,z) = (" 
		  << geom->TPC(t).PlaneLocation(p)[0] << "," 
		  << geom->TPC(t).PlaneLocation(p)[1] << "," 
		  << geom->TPC(t).PlaneLocation(p)[2] << "); \n\t\tpitch from plane 0 is "
		  << geom->TPC(t).Plane0Pitch(p) << "; \n\t\tOrientation "
	          << geom->TPC(t).Plane(p).Orientation() << ", View "
	          << geom->TPC(t).Plane(p).View() << ", Wire angle "
	          << geom->TPC(t).Plane(p).Wire(0).ThetaZ()
		  << std::endl;
      
      std::cout << std::endl << "\t\tTPC Dimensions: " << 2.*geom->TPC(t).HalfWidth()
		<< " x " << 2.*geom->TPC(t).HalfHeight() 
		<< " x " << geom->TPC(t).Length() << std::endl;
      std::cout << std::endl << "\t\tTPC Active Dimensions: " << 2.*geom->TPC(t).ActiveHalfWidth()
		<< " x " << 2.*geom->TPC(t).ActiveHalfHeight() 
		<< " x " << geom->TPC(t).ActiveLength() << std::endl;
      std::cout << "\t\tTPC mass: " << geom->TPC(t).ActiveMass() << std::endl;

      geo::DriftDirection_t dir = geom->TPC(t).DriftDirection();
      if     (dir == geo::kNegX) std::cout << "\t\tdrift direction is towards negative x values" << std::endl;
      else if(dir == geo::kPosX) std::cout << "\t\tdrift direction is towards positive x values" << std::endl;
      else if(dir == geo::kPosY) std::cout << "\t\tdrift direction is towards positive y values" << std::endl;
      else{
	std::cout << "\t\tdrift direction is unknown, assert" << std::endl;
	assert(0);
      }
    }
    
    return;
  }

  //......................................................................
  void GeometryTest::testChannelToWire()
  {
    art::ServiceHandle<geo::Geometry> geom;

    unsigned int t = 0;
    unsigned int p = 0;
    unsigned int w = 0;
    for(unsigned int tpc = 0; tpc < geom->NTPC(); ++tpc){
      for(unsigned int plane = 0; plane < geom->TPC(tpc).Nplanes(); ++plane){
	for(unsigned int wire = 0; wire < geom->TPC(tpc).Plane(plane).Nwires(); ++wire){
	  unsigned int channel = geom->PlaneWireToChannel(plane, wire, tpc);
	  geom->ChannelToWire(channel, t, p, w);

	  if(t != tpc || p != plane || w != wire) assert(0);
	}
      }
    }

    return;
  }

  //......................................................................
  void GeometryTest::testFindPlaneCenters()
  {
    art::ServiceHandle<geo::Geometry> geom;

    double xyz[3] = {0.},   xyzW[3] = {0.};
    for(size_t i = 0; i < geom->Nplanes(); ++i){ 
      geom->Plane(i).LocalToWorld(xyz,xyzW);
      std::cout << "\n\tplane " << i << " is centered at (x,y,z) = (" << xyzW[0] << "," << xyzW[1]
		<< "," << xyzW[2] << ")" << std::endl;
    } 
  } 

  //......................................................................
  void GeometryTest::testWirePos() 
  {
    art::ServiceHandle<geo::Geometry> geom;

    double xyz[3] = {0.};
    for(size_t t = 0; t < geom->NTPC(); ++t){
      const geo::TPCGeo* tpc = &geom->TPC(t); 
      double pitch = tpc->WirePitch();

      for (size_t i=0; i < tpc->Nplanes(); ++i) {
	const geo::PlaneGeo* plane = &tpc->Plane(i);
	double xyzprev[3] = {-1., -1., -1};
	for (size_t j=0; j<plane->Nwires(); ++j) {
	  const geo::WireGeo wire = plane->Wire(j);
	  wire.GetCenter(xyz);
	  if(xyz[2] < xyzprev[2]){
	    std::cout << "\n\twires do not increase in z order" << std::endl;
	    assert(0);
	  }

	  // Test for uniform perpendicular spacing.

	  if(xyzprev[2] > 0. && plane->Orientation() == geo::kVertical) {
	    double dx = xyz[0] - xyzprev[0];
	    if(std::abs(dx) < 1.e-10)
	      dx = 0.;
	    double dy = xyz[1] - xyzprev[1];
	    if(std::abs(dy) < 1.e-10)
	      dy = 0.;
	    double dz = xyz[2] - xyzprev[2];
	    if(std::abs(dz) < 1.e-10)
	      dz = 0.;
	    double s = std::sin(wire.ThetaZ());
	    double c = std::cos(wire.ThetaZ());
	    double du = s * dz - c * dy;
	    if(std::abs(du) < 1.e-10)
	      du = 0.;
	    double dv = c * dz + s * dy;
	    if(std::abs(dv) < 1.e-10)
	      dv = 0.;
	    /*
	    std::cout << "tpc=" << t 
		      << ", plane=" << i
		      << ", wire=" << j 
		      << ", dx=" << dx
		      << ", dy=" << dy
		      << ", dz=" << dz
		      << ", du=" << du
		      << ", dv=" << dv << std::endl;
	    */
	    if(std::abs(du - pitch) > 1.e-6) {
	      std::cout << "\ntpc=" << t 
			<< ", plane=" << i
			<< ", wire=" << j 
			<< ", dx=" << dx
			<< ", dy=" << dy
			<< ", dz=" << dz
			<< ", du=" << du
			<< ", dv=" << dv << std::endl;
	      std::cout << "Wire spacing does not match pitch." << std::endl;
	      assert(0);
	    }
	  }
	  xyzprev[0] = xyz[0];
	  xyzprev[1] = xyz[1];
	  xyzprev[2] = xyz[2];
	}
      }// end loop over planes
    }// end loop over tpcs

    std::cout << "\n\ttesting closest channel algorithm..." << std::endl;
    // Even if you comment it out, please leave the TStopWatch code
    // in this code for additional testing. The NearestChannel routine
    // is the most frequently called in the simulation, so its execution time
    // is an important component of LArSoft's speed.
    TStopwatch stopWatch;
    stopWatch.Start();

    // get a wire and find its center
    for(unsigned int t = 0; t < geom->NTPC(); ++t){
      for(unsigned int p = 0; p < geom->Nplanes(t); ++p){
	for(unsigned int w = 0; w < geom->Plane(p).Nwires(); ++w){
	
	  const geo::WireGeo& wire = geom->TPC(t).Plane(p).Wire(w);
	  const double pos[3] = {0., 0.0, 0.};
	  double posWorld[3] = {0.};
	  wire.LocalToWorld(pos, posWorld);

	  unsigned int t1 = 0; 
	  unsigned int p1 = 0; 
	  unsigned int w1 = 0;
	  unsigned int nearest = 0;
	  try{
	    nearest = geom->NearestChannelFast(posWorld,p,t);
	    //nearest = geom->NearestChannel(posWorld);
	  }
	  catch(cet::exception &e){
	    mf::LogWarning("GeoTestCaughtException") << e;
	  }
	  try{
	    geom->ChannelToWire(nearest, t1, p1, w1);
	  }
	  catch(cet::exception &e){
	    mf::LogWarning("GeoTestCaughtException") << e;
	  }
	  if(t != t1 || p != p1 || w != w1){
	    std::cout << "test point is at " << posWorld[0] << " " << posWorld[1] << " " << posWorld[2] << "\n"
		      << "nearest channel is " << nearest << " for " << t << " " << p << " " << w << "\n"
		      << "tpc/plane/wire is " << t1 << " " << p1 << " " << w1 << std::endl;
	    assert(0);
	  }
	} // end loop over wires
      } // end loop over planes
    }// end loop over tpcs

    stopWatch.Stop();
    std::cout << "\tdone testing closest channel" << std::endl;
    stopWatch.Print();

    // trigger an exception with NearestChannel
    std::cout << "\tattempt to cause an exception to be caught when looking for a nearest channel" << std::endl;

    double posWorld[3] = {geom->CryostatHalfWidth(),
			  geom->CryostatHalfHeight(),
			  geom->CryostatLength()};

    try{
      geom->NearestChannel(posWorld);
    }
    catch(cet::exception &e){
      mf::LogWarning("GeoTestCaughtException") << e;
    }

  }

  //......................................................................
  void GeometryTest::testWirePitch()
  {
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all planes and wires to be sure the pitch is consistent

    // hard code the value we think it should be for each detector
    double shouldbe = 0.4; // true for ArgoNeuT
    if(geom->DetId() == geo::kMicroBooNE) shouldbe = 0.3;
    else if(geom->DetId() == geo::kLBNE)  shouldbe = 0.5;
    else if(geom->DetId() == geo::kBo)    shouldbe = 0.3;

    for(size_t t = 0; t < geom->NTPC(); ++t){
      for(size_t p = 0; p < geom->TPC(t).Nplanes(); ++p){
	for(size_t w = 0; w < geom->TPC(t).Plane(p).Nwires()-1; ++w){
	  double pitch = geom->TPC(t).WirePitch(w, w+1, p);
	  // get the wire pitch
	  if(fabs(pitch - shouldbe) > 0.1*shouldbe){
	    std::cout << "\n\tunexpected pitch: " << pitch << "/" << shouldbe 
		      << " for plane: " << p 
		      << " wires: " << w << ", " << w+1 << std::endl;
	    assert(0);
	  }// end if pitch is wrong
	}// end loop over wires
      }// end loop over planes
    }// end loop over TPCs

  }

  //......................................................................
  void GeometryTest::testPlanePitch()
  {
    art::ServiceHandle<geo::Geometry> geom;

    // loop over all planes to be sure the pitch is consistent

    // hard code the value we think it should be for each detector
    double shouldbe = 0.4; // true for ArgoNeuT
    if(geom->DetId() == geo::kMicroBooNE) shouldbe = 0.3;
    else if(geom->DetId() == geo::kLBNE)  shouldbe = 0.5;
    else if(geom->DetId() == geo::kBo)    shouldbe = 0.3;

    for(size_t t = 0; t < geom->NTPC(); ++t){
      for(size_t p = 0; p < geom->TPC(t).Nplanes()-1; ++p){
	double pitch = fabs(geom->TPC(t).PlanePitch(p, p+1));
	if(fabs(pitch - shouldbe) > 0.1*shouldbe){
	  std::cout << "\n\tunexpected pitch: " << pitch << "/" << shouldbe 
		    << std::endl;
	  assert(0);
	}// end if wrong pitch
      }// end loop over planes
    }// end loop over TPCs

  }

  //......................................................................

  void GeometryTest::testStepping()
  {
    art::ServiceHandle<geo::Geometry> geom;

    //
    // Test stepping. Example is similar to what one would do for photon
    // transport. Rattles photons around inside the scintillator
    // bouncing them off walls.
    //
    double xyz[3],   xyzWire[3] = {0,0,0};
    double dxyz[3], dxyzWire[3] = {0,sin(0.1),cos(0.1)};

    geom->Plane(1).Wire(0).LocalToWorld(xyzWire,xyz);
    geom->Plane(1).Wire(0).LocalToWorldVect(dxyzWire,dxyz);

    std::cout << "\n\t" << xyz[0]  << "\t" << xyz[1]  << "\t" << xyz[2]  << std::endl;
    std::cout << "\t"   << dxyz[0] << "\t" << dxyz[1] << "\t" << dxyz[2] << std::endl;

    gGeoManager->InitTrack(xyz, dxyz);
    for (int i=0; i<10; ++i) {
      const double* pos = gGeoManager->GetCurrentPoint();
      const double* dir = gGeoManager->GetCurrentDirection();
      std::cout << "\tnode = " << 
	gGeoManager->GetCurrentNode()->GetName() << std::endl;
      std::cout << "\t\tpos=" << "\t"
		<< pos[0] << "\t"
		<< pos[1] << "\t"
		<< pos[2] << std::endl;
      std::cout << "\t\tdir=" << "\t"
		<< dir[0] << "\t"
		<< dir[1] << "\t"
		<< dir[2] << std::endl;
      std::cout << "\t\tmat = " << 
	gGeoManager->GetCurrentNode()->
	GetVolume()->GetMaterial()->GetName() << std::endl;
      
      gGeoManager->FindNextBoundary();
      /* const double* norm = */ gGeoManager->FindNormal();
      gGeoManager->Step(kTRUE,kTRUE);
      
    }

    xyz[0] = 306.108; xyz[1] = -7.23775; xyz[2] = 856.757;
    gGeoManager->InitTrack(xyz, dxyz);
    std::cout << "\tnode = " << 
      gGeoManager->GetCurrentNode()->GetName() << std::endl;
    std::cout << "\tmat = " << 
      gGeoManager->GetCurrentNode()->
      GetVolume()->GetMaterial()->GetName() << std::endl;

    gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->Print();

  }

  //......................................................................

  void GeometryTest::testProject() 
  {
    art::ServiceHandle<geo::Geometry> geom;

    double xlo, xhi;
    double ylo, yhi;
    double zlo, zhi;
    geom->WorldBox(&xlo, &xhi, &ylo, &yhi, &zlo, &zhi);
  
    double xyz[3]  = {0.0, 0.0, 0.0};
    double dxyz1[3] = { 1.0, 0.0, 0.0};
    double dxyz2[3] = {-1.0, 0.0, 0.0};
    double dxyz3[3] = { 0.0, 1.0, 0.0};
    double dxyz4[3] = { 0.0,-1.0, 0.0};
    double dxyz5[3] = { 0.0, 0.0, 1.0};
    double dxyz6[3] = { 0.0, 0.0,-1.0};

    double xyzo[3];
    geo::ProjectToBoxEdge(xyz, dxyz1, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (fabs(xyzo[0]-xhi)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz2, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (fabs(xyzo[0]-xlo)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz3, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (fabs(xyzo[1]-yhi)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz4, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (fabs(xyzo[1]-ylo)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz5, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (fabs(xyzo[2]-zhi)>1.E-6) abort();

    geo::ProjectToBoxEdge(xyz, dxyz6, xlo, xhi, ylo, yhi, zlo, zhi, xyzo);
    if (fabs(xyzo[2]-zlo)>1.E-6) abort();
  }


}//end namespace
