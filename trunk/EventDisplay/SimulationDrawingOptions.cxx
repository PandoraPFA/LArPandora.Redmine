////////////////////////////////////////////////////////////////////////
// $Id: SimulationDrawingOption.cxx,v 1.16 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the simulation
//sss
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include "EventDisplay/SimulationDrawingOptions.h"
#include <iostream>

namespace evd {

  //......................................................................
  SimulationDrawingOptions::SimulationDrawingOptions(fhicl::ParameterSet const& pset, 
						     art::ActivityRegistry& reg)
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  SimulationDrawingOptions::~SimulationDrawingOptions() 
  {
  }

  //......................................................................
  void SimulationDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fShowMCTruthText         = pset.get< bool        >("ShowMCTruthText",         true);
    fShowMCTruthVectors      = pset.get< bool        >("ShowMCTruthVectors",      true);
    fShowMCTruthTrajectories = pset.get< bool        >("ShowMCTruthTrajectories", true);
    fMinEnergyDeposition     = pset.get< double      >("MinimumEnergyDeposition"      );
    fG4ModuleLabel           = pset.get< std::string >("G4ModuleLabel"                ); 
  }
  
}
////////////////////////////////////////////////////////////////////////
