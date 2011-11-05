///////////////////////////////////////////////////////////////////////
/// \file    CheckBackTracking.h
/// \brief   test module for confirming that the MC backtracking works
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef CHEAT_CHECKBACKTRACKING_H
#define CHEAT_CHECKBACKTRACKING_H
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"

namespace cheat {
  class CheckBackTracking : public art::EDAnalyzer {
  public:
    explicit CheckBackTracking(fhicl::ParameterSet const& pset);
    virtual ~CheckBackTracking();

    void analyze(art::Event const& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fHitModuleLabel;    ///< label for module creating recob::Hit objects	   
    std::string fG4ModuleLabel;     ///< label for module running G4 and making particles, etc

  };
}
#endif
