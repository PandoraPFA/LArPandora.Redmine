////////////////////////////////////////////////////////////////////////
//
// Display parameters for Argoneut data
//
// \author soderber@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef ARGONEUTDRAWINGOPTIONS_H
#define ARGONEUTDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace argoevd {
  class ArgoneutDrawingOptions 
  {
  public:
    ArgoneutDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~ArgoneutDrawingOptions();

    void reconfigure(fhicl::ParameterSet const& pset);
    
    std::string          fMINOSLabel;               ///< label for MINOS track collection
    std::string    	 fMatchLabel;               ///< label for T962/MINOS track associations			  
   
  };
}//namespace
#endif // __CINT__
#endif

