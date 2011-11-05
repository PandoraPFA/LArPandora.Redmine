////////////////////////////////////////////////////////////////////////
// $Id: RawDrawingOption.cxx,v 1.16 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data 
//
// \author jpaley@indiana.edu
////////////////////////////////////////////////////////////////////////
#include "EventDisplay/RawDrawingOptions.h"
#include <iostream>

namespace evd {

  //......................................................................
  RawDrawingOptions::RawDrawingOptions(fhicl::ParameterSet const& pset, 
				       art::ActivityRegistry& reg) 
  {
    this->reconfigure(pset);
  }
  
  //......................................................................
  RawDrawingOptions::~RawDrawingOptions() 
  {
  }

  void RawDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDrawRawDataOrCalibWires    = pset.get< int         >("DrawRawDataOrCalibWires");
    fScaleDigitsByCharge     	= pset.get< int         >("ScaleDigitsByCharge");    
    fTicksPerPoint              = pset.get< int         >("TicksPerPoint");	  
    fMinSignal                  = pset.get< double      >("MinimumSignal");	  
    fTicks                      = pset.get< double      >("TotalTicks", 2048);	  
    fAxisOrientation         	= pset.get< int         >("AxisOrientation", 0);     
    fRawDataLabel               = pset.get< std::string >("RawDataLabel","daq");
    fTPC                        = pset.get< unsigned int>("TPC", 0);
  }  
}
////////////////////////////////////////////////////////////////////////
