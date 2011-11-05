////////////////////////////////////////////////////////////////////////
// $Id: EvdLayoutOptions.h,v 1.15 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVDLAYOUTOPTIONS_H
#define EVDLAYOUTOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace evd {
  class EvdLayoutOptions 
  {
  public:
    EvdLayoutOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~EvdLayoutOptions();

    void reconfigure(fhicl::ParameterSet const& pset);
    
    int          fShowSideBar;		       ///< 1 to show, 0 don't show
    int    	 fAutoZoomInterest;             ///< Set the automatic zoom to the interest region		  
    int    	 fPrintTotalCharge;             ///< Print out the total charge in an event
    };
}//namespace
#endif // __CINT__
#endif

