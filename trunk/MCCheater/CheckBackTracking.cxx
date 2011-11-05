////////////////////////////////////////////////////////////////////////
/// \file    CheckBackTracking.cxx
/// \brief   check the backtracking of reconstructed hits to simulated particles
/// \version $Id: GeometryTest.cxx,v 1.1 2011/02/17 01:45:48 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <vector>

// ROOT includes

// LArSoft includes
#include "Geometry/geo.h"
#include "MCCheater/BackTracker.h"
#include "MCCheater/CheckBackTracking.h"
#include "RecoBase/recobase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cheat{

  //--------------------------------------------------------------------
  CheckBackTracking::CheckBackTracking(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //--------------------------------------------------------------------
  CheckBackTracking::~CheckBackTracking()
  {
  }

  //--------------------------------------------------------------------
  void CheckBackTracking::reconfigure(fhicl::ParameterSet const& pset)
  {
    fHitModuleLabel    = pset.get< std::string >("HitModuleLabel",    "ffthit"  );
    fG4ModuleLabel     = pset.get< std::string >("G4ModuleLabel",     "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void CheckBackTracking::analyze(art::Event const& evt)
  {

    // grab the sim::ParticleList
    sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);

    // print the list of particles first
    mf::LogInfo("CheckBackTracking") << plist;

    // grab the hits that have been reconstructed
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fHitModuleLabel, hitcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitcol);
    
    // get the sim::SimChannels as well
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);
    
    //now make a vector where each channel in the detector is an 
    // entry
    art::ServiceHandle<geo::Geometry> geo;
    std::vector<const sim::SimChannel*> scs(geo->Nchannels(),0);

    for(size_t i = 0; i < sccol.size(); ++i) scs[sccol[i]->Channel()] = sccol[i];

    // loop over the hits and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();

    // adopt an EmEveIdCalculator to find the eve ID.  
    // will return a primary particle if it doesn't find 
    // a responsible particle for an EM process
    plist.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
    
    // make a collection of the distinct eve ID values
    std::set<int> eveIDs;

    while( itr != hits.end() ){

      // print the truth information for this hit
      mf::LogInfo("CheckBackTracking") << *((*itr).get()) << "\n channel is: " << (*itr)->Channel()
				       << " " << scs[(*itr)->Channel()];

      std::vector<cheat::TrackIDE> trackides = cheat::BackTracker::HitToTrackID(*(scs[(*itr)->Channel()]), *itr);
      std::vector<cheat::TrackIDE> eveides   = cheat::BackTracker::HitToEveID(plist, *(scs[(*itr)->Channel()]), *itr);
      std::vector<double>          xyz       = cheat::BackTracker::HitToXYZ(*(scs[(*itr)->Channel()]), *itr);

      mf::LogInfo("CheckBackTracking") << "hit weighted mean position is (" 
				       << xyz[0] << "," << xyz[1] << "," << xyz[2] << ")";
	
      for(size_t t = 0; t < trackides.size(); ++t){

	// find the Eve particle for the current trackID
	int eveID = plist.EveId( trackides[t].trackID );

	mf::LogInfo("CheckBackTracking") << "track id: " << trackides[t].trackID 
					 << " contributed " << trackides[t].energy << "/" 
					 << trackides[t].energyFrac 
					 << " to the current hit and has eveID: "
					 << eveID;
      }

      for(size_t e = 0; e < eveides.size(); ++e){
	mf::LogInfo("CheckBackTracking") << "eve id: " << eveides[e].trackID 
					 << " contributed " << eveides[e].energy << "/" 
					 << eveides[e].energyFrac 
					 << " to the current hit";

	if(eveIDs.find(eveides[e].trackID) == eveIDs.end()) eveIDs.insert(eveides[e].trackID);
      }

      itr++;
    }// end loop over hits

    // now check the calculation of the purity and efficiency for hit collections
    // first have to make an art::PtrVector<recob::Hit>
    art::PtrVector<recob::Hit> hitPtrVec;
    for(size_t h = 0; h < hits.size(); ++h) hitPtrVec.push_back(hits[h]);

    // loop over the eveID values and calculate the purity and efficiency for each
    std::set<int>::iterator setitr = eveIDs.begin();
    while( setitr != eveIDs.end() ){

      std::set<int> id;
      id.insert(*setitr);
      mf::LogInfo("CheckBackTracking") << "eve ID: " << *setitr 
				       << " purity: " 
				       << cheat::BackTracker::HitCollectionPurity(plist, id, scs, hitPtrVec)
				       << " efficiency: "
				       << cheat::BackTracker::HitCollectionEfficiency(plist, id, scs, hitPtrVec, hitPtrVec, geo::k3D);

      
      setitr++;
    }

    return;

  } // end analyze

} // end namespace
