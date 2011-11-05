////////////////////////////////////////////////////////////////////////
/// \file  Particle.cxx
/// \brief Description of a particle passed to Geant4
///
/// \version $Id: EveIdCalculator.cxx,v 1.1 2010/05/13 16:12:20 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
#include "Simulation/sim.h"
#include "SimulationBase/simbase.h"

#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include <iterator>
#include <iostream>
#include <climits>

namespace sim {

  // static variables

  // How do we indicate an uninitialized variable?  I don't want to
  // use "0" for PDG, because that's potentially a valid value.
  // Instead, let the compiler give us a value.  The following
  // template (from climits) evaluates the lower possible negative
  // number that you can store in an int.

  const int Particle::s_uninitialized = std::numeric_limits<int>::min();
  
  Particle::Particle()
  {
  }

  // Standard constructor.
  Particle::Particle(const int trackId, 
		     const int pdg, 
		     const art::Ptr<simb::MCTruth>& primary,
// 		     const unsigned int primaryIndex,
		     const std::string process,
		     const int mother, 
		     const double mass,
		     const int status)
    : simb::MCParticle(trackId, pdg, process, mother, mass, status)
    , fprimary(primary)
//     , fPrimaryIndex(s_uninitialized)
  {
  }

  Particle::~Particle() 
  {
  }

  // Copy constructor.  Note that since this class inherits from
  // TObject, we have to copy its information explicitly.
  Particle::Particle( const Particle& rhs ) 
  {
    fstatus       = rhs.StatusCode();
    ftrackId      = rhs.TrackId();
    fpdgCode      = rhs.PdgCode();
    fmother       = rhs.Mother();
    fprocess      = rhs.Process();
    ftrajectory   = rhs.Trajectory();
    fmass         = rhs.Mass();
    fpolarization = rhs.Polarization();
    fWeight       = rhs.Weight();

    for(int i = 0; i < rhs.NumberDaughters(); ++i)
      fdaughters.insert(rhs.Daughter(i));
    fprimary      = rhs.fprimary; // This is a art::Ptr<> not a real pointer.
  }

  // Return the position of the last point in the trajectory list.
  const TVector3 Particle::EndPoint() const
  {
    // If the trajectory is empty, return a zero vector.
    if ( ftrajectory.empty() ) return TVector3();
    return (*ftrajectory.rbegin()).first.Vect();
  }

} // namespace sim
