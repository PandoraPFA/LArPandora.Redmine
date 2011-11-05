////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parameters.h
/// \brief Store parameters for running LArG4
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// When a particular derived configuration is specified in a job configuration document, the 
// only place where that information is loaded into LArSoft is via the job control class.
// As an example, say we specify we want to use the "optical" derived configuration for LArG4,
// only the LArG4 object can see that we chose that configuration and hence load the appropriate
// set of parameters.
//
// This class exists to pass parameters loaded by the job controlling object to other 
// related classes. To load or set a parameter, create an instance of LArG4Parameters and
// use the getting and setting methods, all of which point to one particular instance of
// the class, TheLArG4Parameters
//
// Note- I plan to come back to this class andmake it more general (ie load the whole
// config and not use specific get / set methods) soon.
//
// Ben Jones, MIT, March 2010


#include <string>

#ifndef LArG4Parameters_h
#define LArG4Parameters_h 1

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Utilities/PhysicalConstants.h"

namespace sim {

  class LArG4Parameters
    {
    public:
      LArG4Parameters(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~LArG4Parameters() {}
      
      void reconfigure(fhicl::ParameterSet const& pset);

      const int    OpVerbosity()                       const { return fOpVerbosity;            }
      const double ParticleKineticEnergyCut() 	       const { return fParticleKineticECut;    }
      const bool   StoreTrajectories()        	       const { return fStoreTrajectories;      }
      const bool   DrawNeutrals()             	       const { return fDrawNeutrals;           }
      const double VisualizationEnergyCut()   	       const { return fVisualizationEnergyCut; }
      const bool   UseCustomPhysics()         	       const { return fUseCustomPhysics;       }
      const double RecombA()                  	       const { return util::kRecombA;          }
      const double Recombk()                  	       const { return util::kRecombk;          }
      const double GeVToElectrons()           	       const { return util::kGeVToElectrons;   }
      const double LongitudinalDiffusion()    	       const { return fLongitudinalDiffusion;  }
      const double TransverseDiffusion()      	       const { return fTransverseDiffusion;    }
      const double ElectronClusterSize()      	       const { return fElectronClusterSize;    }
      const std::vector<std::string>& EnabledPhysics() const { return fEnabledPhysics;         }
      const int    K0Bias()                            const { return fK0Bias;                 }
      const int    MNXBias()                  	       const { return fXBias; 		       }
      const int    MNXSBias()                 	       const { return fXSBias;		       }
      const bool   KeepEMShowerDaughters()             const { return fKeepEMShowerDaughters;  }
      const bool   DisableWireplanes()                 const { return fDisableWireplanes;      }

    private:
      int                      fOpVerbosity;           ///< Verbosity of optical simulation - soon to be depricated
      double                   fParticleKineticECut;   ///< Minimum energy a particle needs before asking Geant4 
                                                       ///< to track it, GeV
      bool                     fStoreTrajectories;     ///< Whether to store full trajectories for every particle 
                                                       ///< simulated by Geant4
      bool                     fDrawNeutrals;          ///< depricated
      double                   fVisualizationEnergyCut;///< depricated, GeV
      bool                     fUseCustomPhysics;      ///< Whether to use a custom list of physics processes or the default
      double                   fLongitudinalDiffusion; ///< Amount of diffusion in the longitudinal direction, cm^2/ns
      double                   fTransverseDiffusion;   ///< Amount of diffusion in the transverse direction, cm^2/ns
      double                   fElectronClusterSize;   ///< Number of ionization electrons in a given cluster 
                                                       ///< to be simulated in the readout simulation
      std::vector<std::string> fEnabledPhysics;        ///< List of enabled physics processes if using Custom physics
      int                      fK0Bias;                ///< Turns on secondary particle bias for K0, Lambda, 
                                                       ///< neutrons in MuNuclear
      int                      fXSBias;                ///< Turns on cross-section bian in MuNuclear
      int                      fXBias;                 ///< Enhancement factor for cross-section bian in MuNuclear, 
                                                       ///< should be <= 100
      bool                     fKeepEMShowerDaughters; ///< Whether to keep the secondary, tertiary, etc. 
                                                       ///< particles from an EM shower in the output
      bool                     fDisableWireplanes;     ///< Turn of LAr sensitivity and remove charge 
                                                       ///< drift simulation - use for running pure optical sims 
    };
}


#endif
