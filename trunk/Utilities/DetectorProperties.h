////////////////////////////////////////////////////////////////////////
// \file DetectorProperties.h
//
// \brief service to contain information about detector electronics, etc
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_DETECTORPROPERTIES_H
#define UTIL_DETECTORPROPERTIES_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


///General LArSoft Utilities
namespace util{
    class DetectorProperties {
    public:

      DetectorProperties(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~DetectorProperties();
      
      void reconfigure(fhicl::ParameterSet const& p);

      const double SamplingRate()   const { return fSamplingRate;   }
      const int    TriggerOffset()  const { return fTriggerOffset;  }
      const double ElectronsToADC() const { return fElectronsToADC; }
      
    private:

      double fSamplingRate;   ///< in ns
      int    fTriggerOffset;  ///< in # of clock ticks
      double fElectronsToADC; ///< conversion factor for # of ionization electrons to 1 ADC count

    }; // class DetectorProperties
} //namespace utils
#endif // UTIL_DETECTOR_PROPERTIES_H
