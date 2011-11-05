////////////////////////////////////////////////////////////////////////
/// \file  SingleGen.cxx
/// \brief Generator for cosmic-rays
///
/// Module designed to produce a set list of particles for a MC event
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <memory>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools includes
#include "SimulationBase/simbase.h"
#include "EventGeneratorBase/evgenbase.h"

// lar includes
#include "EventGenerator/SingleGen.h"
#include "Geometry/geo.h"
#include "SummaryData/summary.h"

#include "TVector3.h"
#include "TDatabasePDG.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace evgen{

  //____________________________________________________________________________
  SingleGen::SingleGen(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    fSeed = pset.get< unsigned int >("Seed", evgb::GetRandomNumberSeed());

    createEngine( fSeed );

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

  }

  //____________________________________________________________________________
  SingleGen::~SingleGen()
  {
  }

  //____________________________________________________________________________
  void SingleGen::reconfigure(fhicl::ParameterSet const& p)
  {
    fMode         = p.get< int                 >("ParticleSelectionMode");
    fPDG          = p.get< std::vector<int>    >("PDG");
    fP0           = p.get< std::vector<double> >("P0");
    fSigmaP       = p.get< std::vector<double> >("SigmaP");
    fPDist        = p.get< int                 >("PDist");
    fX0           = p.get< std::vector<double> >("X0");
    fY0           = p.get< std::vector<double> >("Y0");
    fZ0           = p.get< std::vector<double> >("Z0");
    fSigmaX       = p.get< std::vector<double> >("SigmaX");
    fSigmaY       = p.get< std::vector<double> >("SigmaY");
    fSigmaZ       = p.get< std::vector<double> >("SigmaZ");
    fPosDist      = p.get< int                 >("PosDist");
    fTheta0XZ     = p.get< std::vector<double> >("Theta0XZ");
    fTheta0YZ     = p.get< std::vector<double> >("Theta0YZ");
    fSigmaThetaXZ = p.get< std::vector<double> >("SigmaThetaXZ");
    fSigmaThetaYZ = p.get< std::vector<double> >("SigmaThetaYZ");
    fAngleDist    = p.get< int                 >("AngleDist");

    // do not put fSeed in reconfigure because we don't want to reset 
    // the seed midstream

    return;
  }

  //____________________________________________________________________________
  void SingleGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;

    geo::DetId_t detid = geo->DetId();

    std::auto_ptr<sumdata::RunData> runcol(new sumdata::RunData(detid));

    run.put(runcol);

    return;
  }

  //____________________________________________________________________________
  void SingleGen::produce(art::Event& evt)
  {

    ///auto_ptr allows ownership to be transferred to the art::Event after the put statement
    std::auto_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    Sample(truth);

    truthcol->push_back(truth);

    evt.put(truthcol);

    return;
  }

  //____________________________________________________________________________
  // Draw the type, momentum and position of a single particle from the 
  // FCIHL description
  void SingleGen::SampleOne(unsigned int i, simb::MCTruth &mct){
    //       std::cout << "picking a " << fPDG[i] << std::endl
    // 		<< fP0[i] << " " << fSigmaP[i] << " (" << fX0[i]
    // 		<< "/" << fSigmaX[i] << ", " << fY0[i] 
    // 		<< "/" << fSigmaY[i] << ", " << fZ0[i] 
    // 		<< "/" << fSigmaZ[i] << ") " << fTheta0XZ[i]
    // 		<< "/" << fSigmaThetaXZ[i] << " " << fTheta0YZ[i]
    // 		<< "/" << fSigmaThetaYZ[i] << std::endl;
    //       std::cout << fPDist << " " << fPosDist << " " << fAngleDist << std::endl;

    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    // Choose momentum
    double p = 0.0;
    double m = 0.0;
    if (fPDist == kGAUS) {
      p = gauss.fire(fP0[i], fSigmaP[i]);
    }
    else {
      p = fP0[i] + fSigmaP[i]*(2.0*flat.fire()-1.0);
    }

    //     std::cout << "get the mass" << std::endl;

    static TDatabasePDG  pdgt;
    TParticlePDG* pdgp = pdgt.GetParticle(fPDG[i]);
    if (pdgp) m = pdgp->Mass();
    
    // Choose position
    TVector3 x;
    if (fPosDist == kGAUS) {
      x[0] = gauss.fire(fX0[i], fSigmaX[i]);;
      x[1] = gauss.fire(fY0[i], fSigmaY[i]);
      x[2] = gauss.fire(fZ0[i], fSigmaZ[i]);
    }
    else {
      x[0] = fX0[i] + fSigmaX[i]*(2.0*flat.fire()-1.0);
      x[1] = fY0[i] + fSigmaY[i]*(2.0*flat.fire()-1.0);
      x[2] = fZ0[i] + fSigmaZ[i]*(2.0*flat.fire()-1.0);
    }
    //       std::cout << "set the position "<<std::endl;

    TLorentzVector pos(x[0], x[1], x[2], 0.0);
    
    //       std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;

    // Choose angles
    double thxz, thyz;
    if (fAngleDist == kGAUS) {
      thxz = gauss.fire(fTheta0XZ[i], fSigmaThetaXZ[i]);
      thyz = gauss.fire(fTheta0YZ[i], fSigmaThetaYZ[i]);
    }
    else {
      thxz = fTheta0XZ[i] + fSigmaThetaXZ[i]*(2.0*flat.fire()-1.0);
      thyz = fTheta0YZ[i] + fSigmaThetaYZ[i]*(2.0*flat.fire()-1.0);
    }

    double thxzrad=thxz*M_PI/180.0;	
    double thyzrad=thyz*M_PI/180.0;

    TLorentzVector pvec(p*cos(thyzrad)*sin(thxzrad),p*sin(thyzrad),p*cos(thxzrad)*cos(thyzrad),sqrt(p*p+m*m));
 
 


    int trackid = -1*(i+1); // set track id to -i as these are all primary particles and have id <= 0
    std::string primary("primary");

    simb::MCParticle part(trackid, fPDG[i], primary);
    part.AddTrajectoryPoint(pos, pvec);

    //       std::cout << "add the particle to the primary" << std::endl;

    mct.Add(part);
  }

  //____________________________________________________________________________
  void SingleGen::Sample(simb::MCTruth &mct) 
  {

//     std::cout << "size of particle vector is " << fPDG.size() << std::endl;


    switch (fMode) {
    case 0: // List generation mode: every event will have one of each
	    // particle species in the fPDG array
      for (unsigned int i=0; i<fPDG.size(); ++i) {
	SampleOne(i,mct);
      }//end loop over particles
      break;
    case 1: // Random selection mode: every event will exactly one particle
            // selected randomly from the fPDG array
      {
	art::ServiceHandle<art::RandomNumberGenerator> rng;
	CLHEP::HepRandomEngine &engine = rng->getEngine();
	CLHEP::RandFlat flat(engine);

	unsigned int i=flat.fireInt(fPDG.size());
	SampleOne(i,mct);
      }
      break;
    default:
      mf::LogWarning("UnrecognizeOption") << "SingleGen does not recognize ParticleSelectionMode "
					  << fMode;
      break;
    } // switch on fMode

    return;
  }

}//end namespace evgen
////////////////////////////////////////////////////////////////////////
