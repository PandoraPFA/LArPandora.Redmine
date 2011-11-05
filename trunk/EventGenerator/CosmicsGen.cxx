////////////////////////////////////////////////////////////////////////
/// \file  CosmicsGen.cxx
/// \brief Generator for cosmic-rays
///
/// Module to produce cosmic ray MC using CRY
///
/// \version $Id: CosmicsGen.cxx,v 1.5 2010/04/23 18:37:36 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <getopt.h>
}
#include <memory>

// ROOT includes
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// larsoft includes
#include "SimulationBase/simbase.h"
#include "EventGenerator/CosmicsGen.h"
#include "EventGeneratorBase/evgenbase.h"
#include "Geometry/geo.h"
#include "SummaryData/summary.h"

namespace evgen{

  //____________________________________________________________________________
  CosmicsGen::CosmicsGen(fhicl::ParameterSet const& pset)
    : fCRYHelp(0)
  {
    // Create a random number engine
    int seed = pset.get< unsigned int >("Seed", evgb::GetRandomNumberSeed());
    createEngine(seed);

    this->reconfigure(pset);
    
    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();    
  }

  //____________________________________________________________________________
  CosmicsGen::~CosmicsGen()
  {
    if(fCRYHelp) delete fCRYHelp;
  }

  //____________________________________________________________________________
  void CosmicsGen::reconfigure(fhicl::ParameterSet const& p)
  {
    if(fCRYHelp){
      delete fCRYHelp; 
      fCRYHelp = 0;
    }

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine();

    fCRYHelp = new evgb::CRYHelper(p, engine);

    return;
  }

  //____________________________________________________________________________
  void CosmicsGen::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    fNhitHisto = tfs->make<TH1F>("fNhitHisto",";hits/20 usec;spills",200,0.0,20000.0);

    fDminHisto0 = tfs->make<TH1F>("fDminHisto0",";d (cm);", 100,0.0,500.E2);
    fDminHisto1 = tfs->make<TH1F>("fDminHisto1",";d (cm);", 100,0.0,500.E2);

    fPhotonAngles   = tfs->make<TH2F>("fPhotonAngles",  ";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fPhotonAnglesLo = tfs->make<TH2F>("fPhotonAnglesLo",";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fPhotonAnglesMi = tfs->make<TH2F>("fPhotonAnglesMi",";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fPhotonAnglesHi = tfs->make<TH2F>("fPhotonAnglesHi",";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fPhotonCosQ     = tfs->make<TH1F>("fPhotonCosQ",    ";cosQ;tracks", 50,-1.0,1.0);
    fPhotonEnergy   = tfs->make<TH1F>("fPhotonEnergy",  ";E (GeV)",     5000,0.0,1000.0);

    fElectronAngles   = tfs->make<TH2F>("fElectronAngles",  ";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fElectronAnglesLo = tfs->make<TH2F>("fElectronAnglesLo",";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fElectronAnglesMi = tfs->make<TH2F>("fElectronAnglesMi",";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fElectronAnglesHi = tfs->make<TH2F>("fElectronAnglesHi",";x;y",36,-180.0,180.0,50,-1.0,1.0);
    fElectronCosQ     = tfs->make<TH1F>("fElectronCosQ",    ";cosQ;tracks",50,-1.0,1.0);
    fElectronEnergy   = tfs->make<TH1F>("fElectronEnergy",  ";E (GeV)",    5000,0.0,1000.0);
  
    fMuonAngles   = tfs->make<TH2F>("fMuonAngles",  ";x;y",        36,-180.0,180.0,50,-1.0,1.0);
    fMuonAnglesLo = tfs->make<TH2F>("fMuonAnglesLo",";x;y",        36,-180.0,180.0,50,-1.0,1.0);
    fMuonAnglesMi = tfs->make<TH2F>("fMuonAnglesMi",";x;y",        36,-180.0,180.0,50,-1.0,1.0);
    fMuonAnglesHi = tfs->make<TH2F>("fMuonAnglesHi",";x;y",        36,-180.0,180.0,50,-1.0,1.0);
    fMuonCosQ     = tfs->make<TH1F>("fMuonCosQ",    ";cosQ;tracks",50,-1.0,1.0);
    fMuonEnergy   = tfs->make<TH1F>("fMuonEnergy",  ";E (GeV)",    5000,0.0,1000.0);

  }

  //____________________________________________________________________________
  void CosmicsGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;

    geo::DetId_t detid = geo->DetId();

    std::auto_ptr<sumdata::RunData> runcol(new sumdata::RunData(detid));

    run.put(runcol);

    return;
  }

  //____________________________________________________________________________
  void CosmicsGen::produce(art::Event& evt)
  {
    std::auto_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kCosmicRay);
    fCRYHelp->Sample(truth, 0);

    // fill some histograms about this event
    art::ServiceHandle<geo::Geometry> geom;
    
    
    // loop over particles in the truth object
    for(int i=0; i<truth.NParticles(); ++i){
      simb::MCParticle particle = truth.GetParticle(i);
      const TLorentzVector& v4 = particle.Position();
      const TLorentzVector& p4 = particle.Momentum();
      double x0[3] = {v4.X(),  v4.Y(),  v4.Z()};
      double dx[3] = {p4.Px(), p4.Py(), p4.Pz()};
      double x1[3] = {0.0,     0.0,     0.5*geom->DetLength()};
      double x2[3];
      double d = geo::ClosestApproach(x1, x0, dx, x2);
      if(d > 0.5*geom->DetLength()) {fDminHisto0->Fill(d); } /// \todo really should see if the particle intersects the detector
      else                         {fDminHisto1->Fill(d); }

      TH1F* hCosQ     = 0;
      TH2F* hAngles   = 0;
      TH2F* hAnglesLo = 0;
      TH2F* hAnglesMi = 0;
      TH2F* hAnglesHi = 0;
      TH1F* hEnergy   = 0;
      if (abs(particle.PdgCode())==13) {
	hCosQ     = fMuonCosQ;
	hAngles   = fMuonAngles;
	hAnglesLo = fMuonAnglesLo;
	hAnglesMi = fMuonAnglesMi;
	hAnglesHi = fMuonAnglesHi;
	hEnergy   = fMuonEnergy;
      }
      else if (abs(particle.PdgCode())==22) {
	hCosQ     = fPhotonCosQ;
	hAngles   = fPhotonAngles;
	hAnglesLo = fPhotonAnglesLo;
	hAnglesMi = fPhotonAnglesMi;
	hAnglesHi = fPhotonAnglesHi;
	hEnergy   = fPhotonEnergy;
      }
      else if (abs(particle.PdgCode())==11) {
	hCosQ     = fElectronCosQ;
	hAngles   = fElectronAngles;
	hAnglesLo = fElectronAnglesLo;
	hAnglesMi = fElectronAnglesMi;
	hAnglesHi = fElectronAnglesHi;
	hEnergy   = fElectronEnergy;
      }
      if (hCosQ!=0) {
	/// \todo really should see if the particle intersects the detector
	if(d < 0.5*geom->DetLength()){
	  double cosq = -p4.Py()/p4.P();
	  double phi  = atan2(p4.Pz(),p4.Px());
	  phi *= 180/M_PI;
	  hCosQ->Fill(cosq);
	  hAngles->Fill(phi,cosq);
	  if      (p4.E()<1.0)  hAnglesLo->Fill(phi,cosq);
	  else if (p4.E()<10.0) hAnglesMi->Fill(phi,cosq);
	  else                  hAnglesHi->Fill(phi,cosq);
	  hEnergy->Fill(p4.E());
	}

      }//end if there is a cos(theta) histogram
    }// loop on particles

    truthcol->push_back(truth);
    evt.put(truthcol);
  
    return;
  }
}
////////////////////////////////////////////////////////////////////////
