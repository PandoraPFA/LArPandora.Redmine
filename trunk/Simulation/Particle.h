////////////////////////////////////////////////////////////////////////
/// \file  Particle.h
/// \brief Particle class
/// \version $Id: ParticleHistory.cxx,v 1.1 2010/04/29 15:38:01 seligman Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This class describes a particle created in the detector Monte
/// Carlo simulation.

#ifndef SIM_PARTICLE_H
#define SIM_PARTICLE_H

#include "art/Persistency/Common/Ptr.h"
#include "SimulationBase/MCParticle.h"

#include <TVector3.h>

#include <string>
#include <iostream>

namespace simb { class MCTruth; }

namespace sim {

  class Particle : public simb::MCParticle {
  public:

    // An indicator for an uninitialized variable (see Particle.cxx).
    static const int s_uninitialized; 

    // Standard constructor.  If the mass is not supplied in the
    // argument, then the PDG mass is used.
    Particle();

    Particle(const int trackId, 
	     const int pdg, 
	     const art::Ptr<simb::MCTruth> &primary,
	     //const unsigned int primaryIndex,
	     const std::string process,
	     const int mother = -1, 
	     const double mass = s_uninitialized,
	     const int status = -1);
    
    Particle( const Particle &rhs);

    // Destructor.
    ~Particle();

    // Return a reference to the MCTruth object if this is a primary
    // particle.  If it's not, return 0.  Note: If the MCTruth
    // objects were not read in for this event, this will always
    // return 0.
    const art::Ptr<simb::MCTruth>& Primary()  const { return fprimary;      }

    // \todo add data member for the primary index in the MCTruth collection
    // \todo once we are using ROOT v5.30+.  Cannot bump the ClassVersion number
    // \todo over 14 before then
    //const unsigned int PrimaryIndex()         const { return fPrimaryIndex; }

    // Return the (x,y,z) of the last point of the particle's trajectory.
    const TVector3  EndPoint()                const;

  private:
    art::Ptr<simb::MCTruth> fprimary;      ///< pointer to primary particle definition
    //unsigned int            fPrimaryIndex; ///< index of the MCTruth object in the event
  };

} // namespace sim


#endif // SIM_PARTICLE_H
