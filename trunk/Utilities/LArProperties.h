////////////////////////////////////////////////////////////////////////
// LArProperties.h
//
// Utility LAr functions
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// msoderbe@syr.edu
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef LARPROPERTIES_H
#define LARPROPERTIES_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


///General LArSoft Utilities
namespace util{
    class LArProperties {
    public:
      LArProperties(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~LArProperties();
      
      double DriftVelocity(double Efield, double Temperature); // in cm/us
      double ElectronLifetime();
      const double   Efield() const { return fEfield; }
      const double   Temperature() const { return fTemperature; }
      double BirksCorrectionAmplitude(double ADC,double calFactor);

    private:
      double                    fEfield;           //kV/cm
      double                    fTemperature;      //kelvin
      double                    fElectronlifetime; //microseconds
    }; // class LArProperties
} //namespace utils
#endif // LARPROPERTIES_H
