////////////////////////////////////////////////////////////////////////
//  LArProperties class
//
//  Implentation of functions to return important LAr properties
//
//  maddalena.antonello@lngs.infn.it
//  ornella.palamara@lngs.infn.it
//  msoderbe@syr.edu
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "Utilities/LArProperties.h"

// ROOT includes
#include "TMath.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
util::LArProperties::LArProperties(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg) :
  fEfield          (pset.get< double >("Efield")          ),
  fTemperature     (pset.get< double >("Temperature")     ),
  fElectronlifetime(pset.get< double >("Electronlifetime"))
{

}

//------------------------------------------------
util::LArProperties::~LArProperties() 
{
 
}

//------------------------------------------------------------------------------------//
double util::LArProperties::DriftVelocity(double Efield, double Temperature){
  
  // Dirft Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  //
  // Efield should have units of kV/cm
  // Temperature should have units of Kelvin

  if(Efield<0.5 || Efield>4.0) 
    mf::LogWarning("LArProperties") << "DriftVelocity Warning! : E-field value of " << Efield 
				    << " kV/cm is outside of range covered by drift"
				    << " velocity parameterization.";


  if(Temperature<87.0 || Temperature>94.0) 
    mf::LogWarning("LArProperties") << "DriftVelocity Warning! : Temperature value of " 
				    << Temperature 
				    << " K is outside of range covered by drift velocity"
				    << " parameterization.";

  double P1 = -0.01481; // K^-1
  double P2 = -0.0075;  // K^-1
  double P3 =  0.141;   // (kV/cm)^-1
  double P4 =  12.4;    // kV/cm
  double P5 =  1.627;   // (kV/cm)^-P6
  double P6 =  0.317;
  double T0 =  90.371;  // K

  double vd = ((P1*(Temperature-T0)+1)
	       *(P3*Efield*TMath::Log(1+P4/Efield) + P5*TMath::Power(Efield,P6))
	       +P2*(Temperature-T0));

  vd /= 10.;

  return vd; // in cm/us
}

//------------------------------------------------------------------------------------//
double util::LArProperties::ElectronLifetime()
{

return fElectronlifetime;

}


//-----------------------------------------------------------------------------------//
/*****************************************************/
	
/// The below function assumes that the user has applied the lifetime correction and effective pitch
/// between the wires (usually after 3D reconstruction). Using with mean wire pitch will not give correct results.
/// parameters: dADCdx - Hit amplitude divided by effective pitch for a given 3D track.
/// calFactor - ADC/fC calibration factor - can be obtained from DetProperties as inverse of electronstoADC

double util::LArProperties::BirksCorrectionAmplitude(double dADCdx,double calFactor){
  /// Correction for charge quenching using parameterization from 
  /// S.Amoruso et al., NIM A 523 (2004) 275
  // copied from CaloArgoItaliano.cxx
  
  const static double eCharge=1.6e-4;  // electron charge in fC
  
  double dQdx = dADCdx/(calFactor*eCharge);
  double dEdx;
  float A3t = 0.800; 
  float K3t = 0.0486; // in KV/cm*(g/cm^2)/MeV
  float rho = 1.388; // LAr density in g/cm^3
  double Wion = 23.6e-6 ; //   23.6 eV = 1e, Wion in MeV/e
  double Efield = 0.485;  // Electric Field in the drift region in KV/cm
  K3t = K3t/rho; // KV/MeV
  dEdx = dQdx/(A3t/Wion-K3t/Efield*dQdx); //MeV/cm 
  return dEdx;
  }



