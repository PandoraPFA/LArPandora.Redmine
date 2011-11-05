////////////////////////////////////////////////////////////////////////
// $Id: EvdLayoutOption.cxx,v 1.16 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data 
//
// \author andrzejs@fnal.gov (based on RawDrawingOptions files)
////////////////////////////////////////////////////////////////////////
#include "EventDisplay/EvdLayoutOptions.h"
#include <iostream>

namespace evd {

  //......................................................................
  EvdLayoutOptions::EvdLayoutOptions(fhicl::ParameterSet const& pset, 
				       art::ActivityRegistry& reg) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  EvdLayoutOptions::~EvdLayoutOptions() 
  {
  }

  void EvdLayoutOptions::reconfigure(fhicl::ParameterSet const& pset)
    {
    fShowSideBar    		= pset.get< int         >("ShowSideBar");
    fAutoZoomInterest     	= pset.get< int         >("AutoZoomInterest");  
    fPrintTotalCharge  		= pset.get< int         >("PrintTotalCharge");  
    }  
}
////////////////////////////////////////////////////////////////////////
