////////////////////////////////////////////////////////////////////////
/// \file    ShowerCheater.cxx
/// \brief   create perfectly reconstructed prongs
/// \version $Id: GeometryTest.cxx,v 1.1 2011/02/17 01:45:48 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <vector>

// ROOT includes

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "ShowerFinder/ShowerCheater.h"
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

namespace shwf{

  //--------------------------------------------------------------------
  ShowerCheater::ShowerCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Shower> >();
  }

  //--------------------------------------------------------------------
  ShowerCheater::~ShowerCheater()
  {
  }

  //--------------------------------------------------------------------
  void ShowerCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedClusterLabel = pset.get< std::string >("CheatedClusterLabel", "cluster" );
    fG4ModuleLabel       = pset.get< std::string >("G4ModuleLabel",       "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void ShowerCheater::produce(art::Event& evt)
  {

    // grab the sim::ParticleList
    sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);

    // get the sim::SimChannels
    // get the sim::SimChannels as well
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fG4ModuleLabel, sccol);
    
    //now make a vector where each channel in the detector is an 
    // entry
    art::ServiceHandle<geo::Geometry> geo;
    std::vector<const sim::SimChannel*> scs(geo->Nchannels(),0);
    for(size_t i = 0; i < sccol.size(); ++i) scs[sccol[i]->Channel()] = sccol[i];

    // grab the clusters that have been reconstructed
    art::Handle< std::vector<recob::Cluster> > clustercol;
    evt.getByLabel(fCheatedClusterLabel, clustercol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Cluster> > clusters;
    art::fill_ptr_vector(clusters, clustercol);
    
    // loop over the clusters and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Cluster> >::iterator itr = clusters.begin();

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map< int, std::vector< art::Ptr<recob::Cluster> > > eveClusterMap;
    std::map< int, std::vector< art::Ptr<recob::Cluster> > >::iterator clusterMapItr = eveClusterMap.begin();

    // loop over all clusters and fill in the map
    while( itr != clusters.end() ){

      // in the ClusterCheater module we set the cluster ID to be 
      // the eve particle track ID*1000 + plane*100 + tpc number.  The
      // floor function on the cluster ID / 1000 will give us
      // the eve track ID
      int eveID = floor((*itr)->ID()/1000.);

      clusterMapItr = eveClusterMap.find(eveID);
	
      // is this id already in the map, if so extend the collection 
      // by one hit, otherwise make a new collection and put it in
      // the map
      if( clusterMapItr != eveClusterMap.end() ){
	  ((*clusterMapItr).second).push_back((*itr));
      }
      else{
	std::vector< art::Ptr<recob::Cluster> > clustervec;
	clustervec.push_back(*itr);
	eveClusterMap[eveID] = clustervec;
      }

      itr++;
    }// end loop over clusters

    // loop over the map and make prongs
    std::auto_ptr< std::vector<recob::Shower> > showercol(new std::vector<recob::Shower>);

    for(clusterMapItr = eveClusterMap.begin(); clusterMapItr != eveClusterMap.end(); clusterMapItr++){

      // separate out the hits for each particle into the different views
      std::vector< art::Ptr<recob::Cluster> > eveClusters( (*clusterMapItr).second );

      art::PtrVector<recob::Cluster> ptrvs;

      std::vector<recob::SpacePoint> spacePoints;

      for(size_t c = 0; c < eveClusters.size(); ++c){
	ptrvs.push_back(eveClusters[c]);
	
	// need to make the space points for this prong
	// loop over the hits for this cluster and make 
	// a space point for each one
	// set the SpacePoint ID to be the cluster ID*10000 
	// + the hit index in the cluster PtrVector of hits
	art::PtrVector<recob::Hit> hits = eveClusters[c]->Hits();
	for(size_t h = 0; h < hits.size(); ++h){
	  art::Ptr<recob::Hit> hit = hits[h];
	  std::vector<double> xyz = cheat::BackTracker::HitToXYZ(*(scs[hit->Channel()]), hit);

	  // make a PtrVector containing just this hit
	  art::PtrVector<recob::Hit> sphit;
	  sphit.push_back(hits[h]);

	  // make the space point and set its ID and XYZ
	  recob::SpacePoint sp(sphit);
	  sp.SetID(eveClusters[c]->ID()*10000 + h);
	  sp.SetXYZ(&xyz[0]);

	  spacePoints.push_back(sp);
	}
      }

      // is this prong electro-magnetic in nature or 
      // hadronic/muonic?  EM --> shower, everything else is a track
      if( abs(plist[(*clusterMapItr).first]->PdgCode()) == 11  ||
	  abs(plist[(*clusterMapItr).first]->PdgCode()) == 22  ||
	  abs(plist[(*clusterMapItr).first]->PdgCode()) == 111 ){

	mf::LogInfo("ShowerCheater") << "prong of " << (*clusterMapItr).first 
				    << " is a shower with pdg code "
				    << plist[(*clusterMapItr).first]->PdgCode();

	// add a prong to the collection.  Make the prong
	// ID be the same as the track ID for the eve particle
	showercol->push_back(recob::Shower(ptrvs, spacePoints));
	showercol->at(showercol->size() - 1).SetID((*clusterMapItr).first);

	// get the direction cosine for the eve ID particle
	// just use the same for both the start and end of the prong
	// \todo Should we make the cheater direction cosines for Prong start and end different?
	const TLorentzVector initmom = plist[(*clusterMapItr).first]->Momentum();
	double dcos[3] = { initmom.Px()/initmom.Mag(),
			   initmom.Py()/initmom.Mag(),
			   initmom.Pz()/initmom.Mag() };

	showercol->back().SetDirection(dcos, dcos);

	mf::LogInfo("ShowerCheater") << "adding shower: \n" 
				     << showercol->back()
				     << "\nto collection.";

      }// end if this is a shower
    } // end loop over the map

    evt.put(showercol);

    return;

  } // end produce

} // end namespace
