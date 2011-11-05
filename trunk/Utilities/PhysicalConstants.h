////////////////////////////////////////////////////////////////////////
/// \file  PhysicalConstants.h
/// \brief Collection of Physical constants used in LArSoft
///
/// 
/// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef UTIL_PHYSICALCONSTANTS_H
#define UTIL_PHYSICALCONSTANTS_H

namespace util{

  // Recombination factor coefficients come from Nucl.Instrum.Meth.A523:275-286,2004
  // R = A/(1 + (dE/dx)*k)
  // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
  // from GeV/voxel width
  // A = 0.800 +/- 0.003
  // k = (0.097+/-0.001) g/(MeVcm^2) 
  static double kRecombA        = 0.8;     ///< see Nucl.Instrum.Meth.A523:275-286,2004
  static double kRecombk        = 0.097;   ///< in g/(MeVcm^{2})

  // Conversion for energy deposited in GeV to number of ionization electrons produced
  static double kGeVToElectrons = 4.237e7; ///< 23.6eV per ion pair, 1e9 eV/GeV

  // More constants
  static double kc    = 2.99792458e10;   ///< cm/s
  static double khbar = 6.58211899e-22;  ///< MeVs

  // Conversion factors
  static double kMeterToCentimeter = 1.e2;                  ///< 1 m = 100 cm
  static double kCentimeterToMeter = 1./kMeterToCentimeter; 
  static double kMeterToKilometer  = 1.e-3;                 ///< 1000 m = 1 km
  static double kKilometerToMeter  = 1./kMeterToKilometer;

  static double keVToMeV           = 1.e-6;                 ///< 1e6 eV = 1 MeV
  static double kMeVToeV           = 1./keVToMeV;
}

#endif //UTIL_PHYSICALCONSTANTS_H
