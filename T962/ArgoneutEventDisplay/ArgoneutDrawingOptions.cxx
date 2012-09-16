////////////////////////////////////////////////////////////////////////
// Display parameters for Argoneut data 
//
// \author soderber@fnal.gov
////////////////////////////////////////////////////////////////////////
#include "T962/ArgoneutEventDisplay/ArgoneutDrawingOptions.h"
#include <iostream>

namespace argoevd {

  //......................................................................
  ArgoneutDrawingOptions::ArgoneutDrawingOptions(fhicl::ParameterSet const& pset, 
                                                 art::ActivityRegistry& reg) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  ArgoneutDrawingOptions::~ArgoneutDrawingOptions() 
  {
  }

  void ArgoneutDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMINOSLabel               = pset.get< std::string >("MINOSLabel","minos");
    fMatchLabel               = pset.get< std::string >("MatchLabel","match");
    fPaddlesLabel             = pset.get< std::string >("PaddlesLabel","paddles");
    fDrawPaddles              = pset.get< bool        >("DrawPaddles",false);
    fCoincidenceTime          = pset.get< int         >("CoincidenceTime",25); 
  }  
}
////////////////////////////////////////////////////////////////////////
