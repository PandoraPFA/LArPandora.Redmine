////////////////////////////////////////////////////////////////////////
//  \file DetectorProperties.cxx 
//
//  \brief Implentation of functions to return important detector properties
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

// LArSoft includes
#include "Utilities/DetectorProperties.h"

namespace util{

  //--------------------------------------------------------------------
  DetectorProperties::DetectorProperties(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
  {
    this->reconfigure(pset);
  }

  //--------------------------------------------------------------------
  DetectorProperties::~DetectorProperties() 
  {
    
  }

  //--------------------------------------------------------------------
  void DetectorProperties::reconfigure(fhicl::ParameterSet const& p)
  {
    fSamplingRate   = p.get< double >("SamplingRate"  );
    fTriggerOffset  = p.get< int    >("TriggerOffset" );
    fElectronsToADC = p.get< double >("ElectronsToADC");

    return;
  }

} // namespace
