////////////////////////////////////////////////////////////////////////
/// \file  LArG4.h
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.h,v 1.11 2010/06/04 21:47:27 bjpjones Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This a module.  It has the following functions:
///
/// - Initialize Geant4 physics, detector geometry, and other
///   processing.
///
/// - Accept sim::MCTruth objects from the MC branch of the FMWK
///   Event structure.
///
/// - Pass the primary particles to the Geant4 simulation to calculate
///   "truth" information for the detector response.
///
/// - Pass the truth information to the DetSim branch of the FMWK event.

#ifndef LARG4_LARG4_H
#define LARG4_LARG4_H 1

#include "art/Framework/Core/EDProducer.h"

#include "G4Base/G4Helper.h"
#include "G4Base/ConvertMCTruthToG4.h"

#include <cstring>

// Forward declarations
class G4RunManager;
class G4UImanager;
class G4VisExecutive;

///Geant4 interface
namespace larg4 {  
 
  // Forward declarations within namespace.
  class LArVoxelListAction;
  class ParticleListAction;
  
  class LArG4 : public art::EDProducer{
  public:
 
    /// Standard constructor and destructor for an FMWK module.
    explicit LArG4(fhicl::ParameterSet const& pset);
    virtual ~LArG4();

    /// The main routine of this module: Fetch the primary particles
    /// from the event, simulate their evolution in the detctor, and
    /// produce the detector response.
    void produce (art::Event& evt); 
    void beginJob();

  protected:
    // These variables are "protected" rather than private, because I
    // can forsee that it may be desirable to derive other simulation
    // routines from this one.

    g4b::G4Helper*       fG4Help;                 ///< G4 interface object					   
    LArVoxelListAction*  flarVoxelListAction;  	  ///< Geant4 user action to accumulate LAr voxel information.
    ParticleListAction*  fparticleListAction;  	  ///< Geant4 user action to particle information.		   

    std::string          fG4PhysListName;         ///< predefined physics list to use if not making a custom one
    std::string          fG4MacroPath;            ///< directory path for Geant4 macro file to be 
                                                  ///< executed before main MC processing.
    bool                 fdumpParticleList;       ///< Whether each event's sim::ParticleList will be displayed.
    bool                 fdumpLArVoxelList;       ///< Whether each event's sim::LArVoxelList will be displayed.
    std::string          fPMTSensitiveVolumeName; ///< Name of volume with which to associate PMT 
                                                  ///< sensitive detectors	
    int                  fSmartStacking;          ///< Whether to instantiate and use class to 
                                                  ///< dictate how tracks are put on stack.	


  private:

    std::string          fGenModuleLabel;         ///< name of the process module label that produced 
                                                  ///< the input simb::MCTruths

  };

} // namespace LArG4

#endif // LARG4_LARG4_H
