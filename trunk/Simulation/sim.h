/// \file sim.h
///
/// \brief Tools and modules for checking out the basics of the Monte Carlo
///
/// \author brebel@fnal.gov
///
/// \version $Id: sim.h,v 1.6 2009/10/05 23:22:58 t962cvs Exp $
///
/// Coding note: Never put '#include "sim.h"' in your code. It would
/// force a dependency in your class on every other class in the
/// Simulation directory; e.g., if I changed the MCTruth code, your
/// code that generates Electrons would re-compile.  This class exists
/// solely as a bookkeeping tool.

#ifndef SIM_H
#define SIM_H

///glom all the includes together
#include "Simulation/LArVoxelData.h"
#include "Simulation/LArVoxelID.h"
#include "Simulation/LArVoxelList.h"
#include "Simulation/Particle.h"
#include "Simulation/ParticleHistory.h"
#include "Simulation/ParticleList.h"
#include "Simulation/SimChannel.h"
#include "Simulation/PMTHit.h"
#include "Simulation/EveIdCalculator.h"
#include "Simulation/EmEveIdCalculator.h"
#include "Simulation/PhotonLibraryParameters.h"

///Monte Carlo Simulation
namespace sim{ 

  unsigned int GetRandomNumberSeed();

  // any track id method returns sim::Particle:NoParticleId, it means the
  // associated particle was too low-energy to be written by the
  // detector Monte Carlo.
  static const int NoParticleId = -999;
  
}

#endif// SIM_H
