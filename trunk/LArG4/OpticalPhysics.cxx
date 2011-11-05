// Adapted for LArSoft by Ben Jones, MIT, Sept 09
//

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#ifndef LARG4_OPTICALPHYSICS_CXX
#define LARG4_OPTICALPHYSICS_CXX 1

#include "LArG4/OpticalPhysics.hh"
#include "LArG4/CustomPhysicsFactory.hh"
#include "LArG4/OpBoundaryProcessSimple.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4AdjointhMultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4VMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4OpticalPhoton.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpWLS.hh"
#include "G4OpRayleigh.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "messagefacility/MessageLogger/MessageLogger.h"

//Register optical physics in custom physics list

namespace larg4 {

  CustomPhysicsFactory<OpticalPhysics> optical_factory("Optical"); 

  //-----------------------------------------------------------
  OpticalPhysics::OpticalPhysics(G4int ver, const G4String& name)
    : G4VPhysicsConstructor(name), verbose(ver)
  {
    G4LossTableManager::Instance();
    mf::LogInfo("OpticalPhysics") << "OBJECT BEING CONSTRUCTED IN OPTICAL PHYSICS";
  }
  
   
  //-----------------------------------------------------------
  OpticalPhysics::~OpticalPhysics()
  {}
  
  //-----------------------------------------------------------  
  void OpticalPhysics::ConstructParticle()
  {
    LOG_DEBUG("OpticalPhysics") << "PARTICLES BEING CONSTRUCTED IN OPTICAL PHYSICS";
    // optical photon
    G4OpticalPhoton::OpticalPhotonDefinition();
    
    // gamma
    G4Gamma::Gamma();
    
    // leptons
    G4Electron::Electron();
    G4Positron::Positron();
    G4MuonPlus::MuonPlus();
    G4MuonMinus::MuonMinus();
    
    // mesons
    G4PionPlus::PionPlusDefinition();
    G4PionMinus::PionMinusDefinition();
    G4KaonPlus::KaonPlusDefinition();
    G4KaonMinus::KaonMinusDefinition();
    
    // barions
    G4Proton::Proton();
    G4AntiProton::AntiProton();
    
    // ions
    G4Deuteron::Deuteron();
    G4Triton::Triton();
    G4He3::He3();
    G4Alpha::Alpha();
    G4GenericIon::GenericIonDefinition();
  }
    
  //-----------------------------------------------------------  
  void OpticalPhysics::ConstructProcess()
  {
    // Add standard EM Processes
    LOG_DEBUG("OpticalPhysics") << "PROCESSES BEING CONSTRUCTED IN OPTICAL PHYSICS";
    
    fTheCerenkovProcess            = new G4Cerenkov("Cerenkov");
    fTheScintillationProcess       = new G4Scintillation("Scintillation");
    fTheAbsorptionProcess          = new G4OpAbsorption();
    fTheRayleighScatteringProcess  = new G4OpRayleigh();
    fTheBoundaryProcess            = new OpBoundaryProcessSimple();
    fTheWLSProcess                 = new G4OpWLS();
    
    
    
    fTheCerenkovProcess->SetMaxNumPhotonsPerStep(700);
    fTheCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    fTheCerenkovProcess->SetTrackSecondariesFirst(false);
    
    fTheScintillationProcess->SetScintillationYieldFactor(1.);
    fTheScintillationProcess->SetTrackSecondariesFirst(false);
    
    // Use Birks Correction in the Scintillation process
    
    G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
    fTheScintillationProcess->AddSaturation(emSaturation);
    
    
    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (fTheCerenkovProcess->IsApplicable(*particle)) {
	pmanager->AddProcess(fTheCerenkovProcess);
	pmanager->SetProcessOrdering(fTheCerenkovProcess,idxPostStep);
	std::cout<<"OpticalPhysics : Cerenkov applicable : " << particleName << std::endl;
      }
      if (fTheScintillationProcess->IsApplicable(*particle)) {
	pmanager->AddProcess(fTheScintillationProcess);
	pmanager->SetProcessOrderingToLast(fTheScintillationProcess, idxAtRest);
	pmanager->SetProcessOrderingToLast(fTheScintillationProcess, idxPostStep);
	std::cout<<"OpticalPhysics : Scintillation applicable : " << particleName << std::endl;
      }
  
     if (particleName == "opticalphoton") {
	G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
	pmanager->AddDiscreteProcess(fTheAbsorptionProcess);
	pmanager->AddDiscreteProcess(fTheRayleighScatteringProcess);
	pmanager->AddDiscreteProcess(fTheBoundaryProcess);
	pmanager->AddDiscreteProcess(fTheWLSProcess);
	std::cout<<"OpticalPhysics : Adding processes to opticalphoton"<<std::endl;
      }
    }
    
  }
}

#endif
