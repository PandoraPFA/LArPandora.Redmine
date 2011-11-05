////////////////////////////////////////////////////////////////////////
/// \file  LArG4.cxx
/// \brief Use Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// LArSoft Includes
#include "LArG4/LArG4.h"
#include "LArG4/LArVoxelReadoutGeometry.h"
#include "LArG4/PhysicsList.h"
#include "LArG4/ParticleListAction.h"
#include "LArG4/PMTSensitiveDetector.h"
#include "LArG4/PMTReadoutGeometry.h"
#include "LArG4/LArStackingAction.h"
#include "LArG4/LArVoxelReadout.h"
#include "Simulation/LArG4Parameters.h"

#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Geometry/geo.h"
#include "G4Base/DetectorConstruction.h"
#include "G4Base/UserActionManager.h"

// G4 Includes
#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4UserRunAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserStackingAction.hh>
#include <G4VisExecutive.hh>
#include <G4VUserPhysicsList.hh>
#include <G4SDManager.hh>
#include <Randomize.hh>
#include <G4SDManager.hh>
#include <G4VSensitiveDetector.hh>
#include <globals.hh>

// ROOT Includes

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

namespace larg4 {

  //----------------------------------------------------------------------
  // Constructor
  LArG4::LArG4(fhicl::ParameterSet const& pset) 
    : fG4Help(0)
    , flarVoxelListAction    (0)
    , fparticleListAction    (0)  
    , fG4PhysListName        (pset.get< std::string >("G4PhysListName","larg4::PhysicsList"))
    , fdumpParticleList      (pset.get< bool        >("DumpParticleList")          )
    , fdumpLArVoxelList      (pset.get< bool        >("DumpLArVoxelList")          )
    , fPMTSensitiveVolumeName(pset.get< std::string >("PMTSensitiveVolumeName")    )
    , fSmartStacking         (pset.get< int         >("SmartStacking")             )	       
    , fGenModuleLabel        (pset.get< std::string >("GenModuleLabel")            )
  {
    LOG_DEBUG("LArG4") << "Debug: LArG4()";

    // get the random number seed, use a random default if not specified    
    // in the configuration file. 
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());
    // setup the random number service for Geant4, the "G4Engine" label is a 
    // special tag setting up a global engine for use by Geant4/CLHEP
    createEngine(seed, "G4Engine");
    
    produces< std::vector<sim::Particle>   >();
    produces< std::vector<sim::PMTHit>     >();
    produces< std::vector<sim::SimChannel> >();

    // constructor decides if initialized value is a path or an environment variable
    cet::search_path sp("FW_SEARCH_PATH");

    sp.find_file(pset.get< std::string >("GeantCommandFile"), fG4MacroPath);
    struct stat sb;
    if (fG4MacroPath.empty() || stat(fG4MacroPath.c_str(), &sb)!=0)
      // failed to resolve the file name
      throw cet::exception("NoG4Macro") << "G4 macro file " 
					<< fG4MacroPath 
					<< " not found!";

  }

  //----------------------------------------------------------------------
  // Destructor
  LArG4::~LArG4() 
  {
    if(fG4Help) delete fG4Help;
  }

  //----------------------------------------------------------------------
  void LArG4::beginJob()
  {
    fG4Help = new g4b::G4Helper(fG4MacroPath,fG4PhysListName);

    // Tell the detector about the parallel LAr voxel geometry.
    std::vector<G4VUserParallelWorld*> pworlds;
    pworlds.push_back( new LArVoxelReadoutGeometry() );
    pworlds.push_back( new PMTReadoutGeometry(fPMTSensitiveVolumeName.c_str()) );

    fG4Help->SetParallelWorlds(pworlds);

    // Use the UserActionManager to handle all the Geant4 user hooks.
    g4b::UserActionManager* uaManager = g4b::UserActionManager::Instance();
    
    // User-action class for accumulating LAr voxels.
    art::ServiceHandle<sim::LArG4Parameters> lgp;
    
    // User-action class for accumulating particles and trajectories
    // produced in the detector.
    fparticleListAction = new ParticleListAction(lgp->ParticleKineticEnergyCut(),
						 lgp->StoreTrajectories(),
						 lgp->KeepEMShowerDaughters());
    uaManager->AddAndAdoptAction(fparticleListAction);

    fG4Help->InitMC();

    // With an enormous detector with lots of rock ala LAr34 (nee LAr20) 
    // we need to be smarter about stacking.

    if (fSmartStacking){
      G4UserStackingAction* stacking_action = new LArStackingAction(fSmartStacking);
      fG4Help->GetRunManager()->SetUserAction(stacking_action);
    }

  }

  //----------------------------------------------------------------------
  void LArG4::produce(art::Event& evt) 
  {
    LOG_DEBUG("LArG4") << "produce()";

    // loop over the lists and put the particles and voxels into the event as collections
    std::auto_ptr< std::vector<sim::Particle>    > partCol  (new std::vector<sim::Particle  >);
    std::auto_ptr< std::vector<sim::SimChannel>  > scCol    (new std::vector<sim::SimChannel>);
    std::auto_ptr< std::vector<sim::PMTHit>      > pmthitCol(new std::vector<sim::PMTHit    >);

    // Fetch the lists of LAr voxels and particles.
    art::ServiceHandle<sim::LArG4Parameters> lgp;

    //look to see if there is any MCTruth information for this
    //event
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    evt.getManyByType(mclists);

    // Need to process Geant4 simulation for each interaction separately.
    // First collect the mc truths into a vector of art::PtrVectors
    std::vector< art::PtrVector<simb::MCTruth> > mcts;
    for(size_t mcl = 0; mcl < mclists.size(); ++mcl){
      
      art::Handle< std::vector<simb::MCTruth> > mclistHandle = mclists[mcl];
      art::PtrVector<simb::MCTruth> truths;
      for(size_t i = 0; i < mclistHandle->size(); ++i){
	art::Ptr<simb::MCTruth> mct(mclistHandle,i);
	truths.push_back(mct);
      }
      mcts.push_back(truths);
    }

    // get the electrons from the LArVoxelReadout sensitive detector
    // Get the sensitive-detector manager.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    
    // Find the sensitive detector with the name "LArVoxelSD".
    G4VSensitiveDetector* sd = sdManager->FindSensitiveDetector("LArVoxelSD");
    PMTSensitiveDetector *thePMTDet = dynamic_cast<PMTSensitiveDetector*>(sdManager->FindSensitiveDetector("PMTSensitiveDetector"));
    
    // If this didn't work, then a sensitive detector with
    // the name "LArVoxelSD" does not exist.
    if ( !sd ){
      mf::LogError("LArG4") << "Sensitive detector 'LArVoxelSD' does not exist";
      assert(0);
    }
    
    // Convert the G4VSensitiveDetector* to a LArVoxelReadout*.
    LArVoxelReadout *larVoxelReadout = dynamic_cast<LArVoxelReadout*>( sd );
    
    // If this didn't work, there is a "LArVoxelSD" detector, but
    // it's not a LArVoxelReadout object.
    if ( !larVoxelReadout ){
      mf::LogError("LArG4") << "Sensitive detector 'LArVoxelSD' is not a "
			    << "LArVoxelReadout object";
      assert(0);
    }

    // now loop over the collection of art::PtrVectors
    for(size_t mc = 0; mc < mcts.size(); ++mc) {

      // The following tells Geant4 to track the particles in this event.
      fG4Help->G4Run(mcts[mc]);

      for(size_t e = 0; e < larVoxelReadout->GetSimChannels().size(); ++e)
	scCol->push_back(larVoxelReadout->GetSimChannels().at(e));

      const sim::ParticleList& particleList = *( fparticleListAction->GetList() );
      for( sim::ParticleList::const_iterator i = particleList.begin(); i != particleList.end(); ++i){
	// copy the particle so that it isnt const
	sim::Particle p(*(*i).second);
	partCol->push_back(p);
      }

      LOG_DEBUG("LArG4") << "now put the SimChannels and Particles in the event";
      
      // Has the user request a detailed dump of the output objects?
      if (fdumpParticleList){
	mf::LogInfo("LArG4") << "Dump sim::ParticleList; size()=" 
			     << particleList.size() << std::endl
			     << particleList;
      }
    
      // Get PMT Hit collection(s) to store from the sensitive detectors and
      // store in the event
      if(thePMTDet){
	const sim::PMTHitCollection& pmtHitCollection = *thePMTDet->GetPMTHitCollection();
	LOG_DEBUG("Optical") << "Storing PMT Hit Collection in Event"; 
	
	for( sim::PMTHitCollection::const_iterator i = pmtHitCollection.begin(); i != pmtHitCollection.end(); ++i){
	  sim::PMTHit ph(*((*i).second));
	  pmthitCol->push_back(ph);
	}
      }// end if there is a PMT detector volume

    }// end loop over interactions

    evt.put(scCol);
    evt.put(partCol);
    evt.put(pmthitCol);

    return;
  }

} // namespace LArG4
