////////////////////////////////////////////////////////////////////////
/// \file  BackTracker.cxx
/// \brief back track the reconstruction to the simulation
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <map>
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes

#include "MCCheater/BackTracker.h"
#include "Utilities/DetectorProperties.h"

namespace cheat{

  //----------------------------------------------------------------------
  BackTracker::BackTracker()
  {
  }
  
  //----------------------------------------------------------------------
  BackTracker::~BackTracker()
  {
  }
  
  //----------------------------------------------------------------------
  const std::vector<TrackIDE> BackTracker::HitToTrackID(sim::SimChannel      const& channel,
							art::Ptr<recob::Hit> const& hit)
  {
    std::vector<TrackIDE> trackIDEs;
    
    double start = hit->StartTime();
    double end   = hit->EndTime();
	
    cheat::BackTracker::ChannelToTrackID(trackIDEs, channel, start, end);

    return trackIDEs;
  }

  //----------------------------------------------------------------------
  // plist is assumed to have adopted the appropriate EveIdCalculator prior to 
  // having been passed to this method. It is likely that the EmEveIdCalculator is
  // the one you always want to use
  const std::vector<TrackIDE> BackTracker::HitToEveID(sim::ParticleList    const& plist,
						      sim::SimChannel      const& channel,
						      art::Ptr<recob::Hit> const& hit)
  {
    std::vector<TrackIDE> eveides;
    std::vector<TrackIDE> trackides = cheat::BackTracker::HitToTrackID(channel, hit);

    // make a map of evd ID values and fraction of energy represented by
    // that eve id in this hit
    std::map<int, float> eveToE;
    
    double totalE = 0.;
    for(size_t t = 0; t < trackides.size(); ++t){
      eveToE[plist.EveId( trackides[t].trackID )] += trackides[t].energy;
      totalE += trackides[t].energy;
    }
    
    // now fill the eveides vector from the map
    for(std::map<int, float>::iterator itr = eveToE.begin(); itr != eveToE.end(); itr++){
      TrackIDE temp;
      temp.trackID    = (*itr).first;
      temp.energyFrac = (*itr).second/totalE;
      temp.energy     = (*itr).second;
      eveides.push_back(temp);
    }

    return eveides;
  }

  //----------------------------------------------------------------------
  // plist is assumed to have adopted the appropriate EveIdCalculator prior to 
  // having been passed to this method. It is likely that the EmEveIdCalculator is
  // the one you always want to use
  std::set<int> BackTracker::GetSetOfEveIDs(sim::ParticleList  const& plist)
  {
    std::set<int> eveIDs;

    sim::ParticleList::const_iterator plitr = plist.begin();
    while(plitr != plist.end() ){
      int eveID = plist.EveId((*plitr).first);
      // look to see if this eveID is already in the set
      if( eveIDs.find(eveID) == eveIDs.end() ) eveIDs.insert(eveID);
      plitr++;
    }

    return eveIDs;
  }

  //----------------------------------------------------------------------
  double BackTracker::HitCollectionPurity(sim::ParticleList                   const& plist,
					  std::set<int>                              trackIDs, 
					  std::vector<const sim::SimChannel*> const& channels,
					  art::PtrVector<recob::Hit>          const& hits)
  {
    // get the list of EveIDs that correspond to the hits in this collection
    // if the EveID shows up in the input list of trackIDs, then it counts
    float total   = 1.*hits.size();;
    float desired = 0.;
    for(size_t h = 0; h < hits.size(); ++h){

      art::Ptr<recob::Hit> hit = hits[h];
      std::vector<TrackIDE> eveIDs = cheat::BackTracker::HitToEveID(plist, *(channels[hit->Channel()]), hit);

      // don't double count if this hit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < eveIDs.size(); ++e){
	if(trackIDs.find(eveIDs[e].trackID) != trackIDs.end()){
	  desired += 1.;
	  break;
	}
      }

    }// end loop over hits

    double purity = 0.;
    if(total > 0) purity = desired/total;

    return purity;
  }

  //----------------------------------------------------------------------
  double BackTracker::HitCollectionEfficiency(sim::ParticleList                   const& plist,
					      std::set<int>                              trackIDs, 
					      std::vector<const sim::SimChannel*> const& channels,
					      art::PtrVector<recob::Hit>          const& hits,
					      art::PtrVector<recob::Hit>          const& allhits,
					      geo::View_t                         const& view)
  {
    // get the list of EveIDs that correspond to the hits in this collection
    // and the energy associated with the desired trackID
    float desired = 0.;
    float total   = 0.;

    for(size_t h = 0; h < hits.size(); ++h){

      art::Ptr<recob::Hit> hit = hits[h];
      std::vector<TrackIDE> eveIDs = cheat::BackTracker::HitToEveID(plist, *(channels[hit->Channel()]), hit);

      // don't double count if this hit has more than one of the
      // desired track IDs associated with it
      for(size_t e = 0; e < eveIDs.size(); ++e){
	if(trackIDs.find(eveIDs[e].trackID) != trackIDs.end()){
	  desired += 1.;
	  break;
	}
      }
    }// end loop over hits

    // now figure out how many hits in the whole collection are associated with this id
    for(size_t h = 0; h < allhits.size(); ++h){

      art::Ptr<recob::Hit> hit = allhits[h];
      std::vector<TrackIDE> eveIDs = cheat::BackTracker::HitToEveID(plist, *(channels[hit->Channel()]), hit);

      // check that we are looking at the appropriate view here
      // in the case of 3D objects we take all hits
      if(hit->View() != view && view != geo::k3D ) continue;

      for(size_t e = 0; e < eveIDs.size(); ++e){
	// don't worry about hits where the energy fraction for the chosen
	// eveID is < 0.1
	// also don't double count if this hit has more than one of the
	// desired track IDs associated with it
	if(trackIDs.find(eveIDs[e].trackID) != trackIDs.end()
	   && eveIDs[e].energyFrac >= 0.1){
	  total += 1.;
	  break;
	}
      }

    }// end loop over all hits
    
    double efficiency = 0.;
    if(total > 0.) efficiency = desired/total;

    return efficiency;
  }

  //----------------------------------------------------------------------
  void BackTracker::ChannelToTrackID(std::vector<TrackIDE>& trackIDEs,
				     sim::SimChannel const& channel,
				     double                 startTime,
				     double                 endTime)
  {
    trackIDEs.clear();

    double totalE = 0.;

    // loop over the electrons in the channel and grab those that are in time 
    // with the identified hit start and stop times
    std::vector<sim::IDE> simides = channel.TrackIDsAndEnergies((unsigned int)startTime, 
								(unsigned int)endTime);

    // first get the total energy represented by all track ids for 
    // this channel and range of tdc values
    for(size_t e = 0; e < simides.size(); ++e)
      totalE += simides[e].numElectrons;
      
      
    // protect against a divide by zero below
    if(totalE < 1.e-5) totalE = 1.;

    // loop over the entries in the map and fill the input vectors
    
    for(size_t e = 0; e < simides.size(); ++e){

      if(simides[e].trackID == sim::NoParticleId) continue;
      
      TrackIDE info;
      info.trackID    = simides[e].trackID;
      info.energyFrac = simides[e].numElectrons/totalE;
      info.energy     = simides[e].numElectrons;

      trackIDEs.push_back(info);

    }

    return;
  }

  //----------------------------------------------------------------------
  // this method assumes that you are passing the appropriate sim::SimChannel
  // in for the specified hit.  Easiest thing to do is make a 
  // std::vector<const sim::SimChannel*>
  // that has a null entry for each channel in the detector 
  // at the start of your module and then fill the 
  // entries that correspond to actual saved sim::SimChannels
  void BackTracker::HitToSimIDEs(sim::SimChannel      const& channel,
				 art::Ptr<recob::Hit> const& hit,
				 std::vector<sim::IDE>&      ides)
  {

    double startTime = hit->StartTime();
    double endTime   = hit->EndTime();

    ides = channel.TrackIDsAndEnergies( (unsigned int) startTime,
					(unsigned int) endTime);
  }

  //----------------------------------------------------------------------
  std::vector<double> BackTracker::SimIDEsToXYZ(std::vector<sim::IDE> const& ides)
  {
    std::vector<double> xyz(3, -999.);

    double x = 0.;
    double y = 0.;
    double z = 0.;
    double w = 0.;

    // loop over electrons.

    for(std::vector<sim::IDE>::const_iterator ie = ides.begin(); ie != ides.end(); ++ie) {

      const sim::IDE& ide = *ie;
      
      double weight = ide.numElectrons;
      
      w += weight;
      x += weight * ide.x;
      y += weight * ide.y;
      z += weight * ide.z;

    }// end loop over sim::IDEs
	
    // if the sum of the weights is still 0, then return
    // the obviously stupid default values
    if(w < 1.e-5) return xyz;

    xyz[0] = x/w;
    xyz[1] = y/w;
    xyz[2] = z/w;

    return xyz;
  }

  //----------------------------------------------------------------------
  std::vector<double> BackTracker::HitToXYZ(sim::SimChannel      const& channel,
					    art::Ptr<recob::Hit> const& hit)
  {
    std::vector<sim::IDE> ides;
    HitToSimIDEs(channel, hit, ides);
    return SimIDEsToXYZ(ides);
  }

  //----------------------------------------------------------------------
  // channels is assumed to have an entry for every channel in the detector,
  // those without actual signals from the MC are assumed to be null pointers
  std::vector<double> BackTracker::SpacePointToXYZ(std::vector<const sim::SimChannel*> const& channels,
						   recob::SpacePoint const& spt)
  {
    // Result vector.

    std::vector<double> xyz(3, 0.);

    // Get hits that make up this space point.

    const art::PtrVector<recob::Hit>& hits = spt.Hits(geo::kU, true);

    // Loop over hits and sum positions of each hit.

    int nhits = hits.size();
    if(nhits > 0) {
      for(art::PtrVector<recob::Hit>::const_iterator i = hits.begin();
	  i != hits.end(); ++i) {
	const art::Ptr<recob::Hit>& phit = *i;

	std::vector<double> xyzhit = HitToXYZ(*(channels[phit->Channel()]), phit);
	xyz[0] += xyzhit[0];
	xyz[1] += xyzhit[1];
	xyz[2] += xyzhit[2];
      }
      xyz[0] /= nhits;
      xyz[1] /= nhits;
      xyz[2] /= nhits;
    }

    // Done.

    return xyz;
  }

} // namespace
