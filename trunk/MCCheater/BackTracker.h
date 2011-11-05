////////////////////////////////////////////////////////////////////////
/// \file  BackTracker.h
/// \brief back track the reconstruction to the simulation
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////


#ifndef CHEAT_BACKTRACKER_H
#define CHEAT_BACKTRACKER_H

#include <vector>
#include "Geometry/geo.h"
#include "Simulation/sim.h"
#include "RawData/raw.h"
#include "RecoBase/recobase.h"

///code to link reconstructed objects back to the MC truth information
namespace cheat{

  typedef struct{
    int trackID;      ///< Geant4 supplied trackID
    float energyFrac; ///< fraction of hit energy from the particle with this trackID
    float energy;     ///< energy from the particle with this trackID
  } TrackIDE;

  class BackTracker
  {

  public:

    BackTracker();
    ~BackTracker();

    // all methods that take a single sim::SimChannel const& assume that
    // the channel for the appropriate hit has already been selected. 
    // The easiest thing to do is make a 
    // std::vector<const sim::SimChannel*>
    // that has a null entry for each channel in the detector at the start
    // of your module and then fill the entries that correspond to actual 
    // saved sim::SimChannels 

    // this method will return the Geant4 track IDs of 
    // the particles contributing ionization electrons to the identified hit
    static const std::vector<TrackIDE> HitToTrackID(sim::SimChannel      const& channel,
						    art::Ptr<recob::Hit> const& hit);
    
    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    static const std::vector<TrackIDE> HitToEveID(sim::ParticleList    const& plist,
						  sim::SimChannel      const& channel,
						  art::Ptr<recob::Hit> const& hit);
    
    // method to return sim::IDE objects associated with a given hit
    static void HitToSimIDEs(sim::SimChannel      const& channel,
			     art::Ptr<recob::Hit> const& hit,
			     std::vector<sim::IDE>&      ides);

    // method to return the XYZ position of the weighted average energy deposition for a given hit
    static std::vector<double> SimIDEsToXYZ(std::vector<sim::IDE> const& ides);

    // method to return the XYZ position of the weighted average energy deposition for a given hit
    static std::vector<double>  HitToXYZ(sim::SimChannel      const& channel,
					 art::Ptr<recob::Hit> const& hit);

    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    static std::vector<double> SpacePointToXYZ(std::vector<const sim::SimChannel*> const& channels,
					       recob::SpacePoint                   const& spt);
    
    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids 
    static double              HitCollectionPurity(sim::ParticleList                   const& plist,
						   std::set<int>                              trackIDs, 
						   std::vector<const sim::SimChannel*> const& channels,
						   art::PtrVector<recob::Hit>          const& hits);
    
    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are 
    // represented in a collection of hits
    static double              HitCollectionEfficiency(sim::ParticleList                   const& plist,
						       std::set<int>                              trackIDs, 
						       std::vector<const sim::SimChannel*> const& channels,
						       art::PtrVector<recob::Hit>          const& hits,
						       art::PtrVector<recob::Hit>          const& allhits,
						       geo::View_t                         const& view);
    
    // method to return all EveIDs corresponding to the given sim::ParticleList
    static std::set<int>       GetSetOfEveIDs(sim::ParticleList const& plist);

  private:

    static void ChannelToTrackID(std::vector<TrackIDE>& trackIDEs,
				 sim::SimChannel const& channel,
				 double                 startTime,
				 double                 endTime);

  };
} // namespace

#endif // CHEAT_BACKTRACK_H
