////////////////////////////////////////////////////////////////////////
/// \file PMTSensitiveDetector.cxx
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Implementation of the PMTSensitiveDetector
//
// See comments in PMTSensitiveDetector.h
//
// Ben Jones, MIT, 06/04/2010
//



#include "LArG4/PMTSensitiveDetector.h"
#include "G4SDManager.hh"

namespace larg4{


  PMTSensitiveDetector::PMTSensitiveDetector(G4String DetectorUniqueName)
    : G4VSensitiveDetector(DetectorUniqueName)
  {
    // Register self with sensitive detector manager
    G4SDManager::GetSDMpointer()->AddNewDetector(this);

    // Get instances of singleton classes
    fThePMTLookup        = PMTLookup::Instance();

    //Create a new hit collection for this sensitive detector and name it
    fThePMTHitCollection = new sim::PMTHitCollection;
    fThePMTHitCollection->SetSDName(DetectorUniqueName);
  }


  sim::PMTHitCollection * PMTSensitiveDetector::GetPMTHitCollection()
  {
    return fThePMTHitCollection;
  }


  G4bool PMTSensitiveDetector::ProcessHits(G4Step * aStep, G4TouchableHistory *)
  {
    sim::PMTPhoton ThePhoton;


    // Get photon data to store in the hit

    ThePhoton.SetInSD      = true;

    ThePhoton.Position     = TLorentzVector(
					    aStep->GetTrack()->GetPosition().x(),
					    aStep->GetTrack()->GetPosition().y(),
					    aStep->GetTrack()->GetPosition().z(),
					    aStep->GetTrack()->GetGlobalTime()
					    );

    ThePhoton.Momentum     = TLorentzVector(
					    aStep->GetTrack()->GetMomentum().x(),
					    aStep->GetTrack()->GetMomentum().y(),
					    aStep->GetTrack()->GetMomentum().z(),
					    aStep->GetTrack()->GetTotalEnergy()
					    );

    // Lookup which PMT we are in
    int PMTID = fThePMTLookup->GetID(aStep->GetPreStepPoint()->GetPhysicalVolume());

    // Add this photon to the relevant PMT hit
    (  *(fThePMTHitCollection->GetHit(PMTID))  ).push_back(ThePhoton);
    (  *(fThePMTHitCollection->GetHit(PMTID))  ).SetID(PMTID);
    
    // Kill this photon track
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    std::cout<<"Particle stepping in PMT" << PMTID <<std::endl;

    return true;
    

  }


  void PMTSensitiveDetector::Initialize(G4HCofThisEvent *)
  {
    fThePMTHitCollection->clear();

  }

}
