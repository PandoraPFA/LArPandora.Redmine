////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parameters.cxx
/// \brief Store parameters for running LArG4
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// Whena particular derived configuration is specifiedin a job configuration document, the
// only place where that information isloaded into LArSoft is via the job control class.
// As an example, say we specify we want to usethe "optical" derived configuration forLArG4,
// only the LArG4 object can see that we chose that configuration and hence load the appropriate
// set of parameters.
//
// This class exists to pass parameters loaded by the job controlling object to other
// related classes. To load or set a parameter, create an instance of LArG4Parameters and
// use the getting and setting methods,all of which point to one particular instance of
// the class, TheLArG4Parameters
//
// Note - I plan to come back to this class and make it more general (ie load the whole
// config and not use specific get / set methods) soon.
//
// Ben Jones, MIT, March 2010



#include "Simulation/LArG4Parameters.h"


namespace sim {

  LArG4Parameters::LArG4Parameters(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    this->reconfigure(pset);
  }

  //--------------------------------------------------------------------------
  void LArG4Parameters::reconfigure(fhicl::ParameterSet const& pset)
  {

    fOpVerbosity            = pset.get< int                      >("OpticalSimVerbosity");     
    fParticleKineticECut    = pset.get< double 			 >("ParticleKineticEnergyCut");
    fStoreTrajectories      = pset.get< bool   			 >("StoreTrajectories");	     
    fDrawNeutrals           = pset.get< bool   			 >("VisualizeNeutrals");	     
    fVisualizationEnergyCut = pset.get< double 			 >("VisualizationEnergyCut");  
    fUseCustomPhysics       = pset.get< bool   			 >("UseCustomPhysics");	     
    fKeepEMShowerDaughters  = pset.get< bool                     >("KeepEMShowerDaughters");
    fLongitudinalDiffusion  = pset.get< double 			 >("LongitudinalDiffusion");   
    fTransverseDiffusion    = pset.get< double 			 >("TransverseDiffusion");     
    fElectronClusterSize    = pset.get< double 			 >("ElectronClusterSize");     
    fEnabledPhysics         = pset.get< std::vector<std::string> >("EnabledPhysics");
    fK0Bias                 = pset.get< int                      >("CosmogenicK0Bias");	     
    fXBias                  = pset.get< int    			 >("CosmogenicXSMNBiasOn");    
    fXSBias                 = pset.get< int     		 >("CosmogenicXSMNBiasFactor");
    fDisableWireplanes      = pset.get< bool                     >("DisableWireplanes");

    // First of last 3 flags above turns on secondary particle bias for 
    // K0s,Lambdas,neutrons in MuNuclear. 
    // The second turns on cross-section bias in MuNuclear.
    // The 3rd is the enhancement factor for XS bias in MuNuclear. Keep it 
    // <=100.

    return;
  }

}
