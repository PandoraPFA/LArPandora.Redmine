///////////////////////////////////////////////////////////////////////
/// \file    GeometryTest.h
/// \brief   unit tests for the Geometry service
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOMETRYTEST_H
#define GEO_GEOMETRYTEST_H
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"

namespace geo { class Geometry; }

///tracking algorithms
namespace geo {
  class GeometryTest : public art::EDAnalyzer {
  public:
    explicit GeometryTest(fhicl::ParameterSet const& pset);
    virtual ~GeometryTest();

    void analyze(art::Event const& evt);
    
  private:

    void testTPC();
    void testChannelToWire();
    void testFindPlaneCenters();
    void testProject();
    void testWirePitch();
    void testPlanePitch();
    void testWirePos();
    void testStepping();
  };
}
#endif
