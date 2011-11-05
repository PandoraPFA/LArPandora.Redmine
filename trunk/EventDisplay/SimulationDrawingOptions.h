////////////////////////////////////////////////////////////////////////
// $Id: SimulationDrawingOption.h,v 1.15 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SIMULATIONDRAWINGOPTIONS_H
#define SIMULATIONDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd {
  class SimulationDrawingOptions 
  {
  public:
    SimulationDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~SimulationDrawingOptions();
    
    void reconfigure(fhicl::ParameterSet const& pset);

    bool        fShowMCTruthText;
    bool        fShowMCTruthVectors;
    bool        fShowMCTruthTrajectories;
    double      fMinEnergyDeposition;
    std::string fG4ModuleLabel;           ///< module label producing sim::SimChannel objects

  };
}//namespace
#endif // __CINT__
#endif

