////////////////////////////////////////////////////////////////////////
// Class:       RecoCheckAna
// Module Type: analyzer
// File:        RecoCheckAna.h
//
// Generated at Fri Jul 15 09:54:26 2011 by Brian Rebel using artmod
// from art v0_07_04.

#include "TH1.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "MCCheater/RecoCheckAna.h"
#include "MCCheater/BackTracker.h"
#include "Simulation/SimListUtils.h"

//-------------------------------------------------------------------
cheat::RecoCheckAna::RecoCheckAna(fhicl::ParameterSet const &p)
{
  this->reconfigure(p);
}

//-------------------------------------------------------------------
cheat::RecoCheckAna::~RecoCheckAna() 
{
  // Clean up dynamic memory and other resources here.
}

//-------------------------------------------------------------------
// We want to determine how well an event finding algorithm works by 
// determining the following:
// 1. How much of the actual energy for a given interaction is 
//    found in the reconstructed event
// 2. How many of the particles from a given interaction are assigned 
//    to the correct reconstructed event
// 3. How close is the reconstructed event primary vertex to the 
//    simulated interaction primary vertex
// 4. What fraction of energy in the reconstructed event is from 
//    particles not associated with the simulated interaction

void cheat::RecoCheckAna::analyze(art::Event const &e) 
{

  // check that this is MC, stop if it isn't
  if(e.isRealData()){
    mf::LogWarning("RecoVetter") << "attempting to run MC truth check on "
				 << "real data, bail";
    return;
  }

  // grab the sim::ParticleList
  sim::ParticleList plist = sim::SimListUtils::GetParticleList(e, fG4ModuleLabel);

  // adopt an EmEveIdCalculator to find the eve ID.  
  // will return a primary particle if it doesn't find 
  // a responsible particle for an EM process
  plist.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
  
  // get the sim::SimChannels 
  std::vector<const sim::SimChannel*> sccol;
  e.getView(fG4ModuleLabel, sccol);
    
  // now make a vector where each channel in the detector is an 
  // entry
  art::ServiceHandle<geo::Geometry> geo;
  std::vector<const sim::SimChannel*> scs(geo->Nchannels(),0);
  for(size_t i = 0; i < sccol.size(); ++i) scs[sccol[i]->Channel()] = sccol[i];

  // get all hits in the event to figure out how many there are
  art::Handle< std::vector<recob::Hit> > hithdl;
  e.getByLabel(fHitModuleLabel, hithdl);
  art::PtrVector<recob::Hit> allhits;
  for(size_t h = 0; h < hithdl->size(); ++h){
    art::Ptr<recob::Hit> hit(hithdl, h);
    allhits.push_back(hit);
  }

  // define variables to hold the reconstructed objects
  art::Handle< std::vector<recob::Cluster> > clscol;
  std::vector< art::Ptr<recob::Cluster> >    clusters;
  art::Handle< std::vector<recob::Track> >   trkcol;
  std::vector< art::Ptr<recob::Prong> >      tracks;
  art::Handle< std::vector<recob::Shower> >  shwcol;
  std::vector< art::Ptr<recob::Prong> >      showers;
  art::Handle< std::vector<recob::Vertex> >  vtxcol;
  std::vector< art::Ptr<recob::Vertex> >     vertices;
  art::Handle< std::vector<recob::Event> >   evtcol;
  std::vector< art::Ptr<recob::Event> >      events;

  if(fCheckClusters){
    try{
      e.getByLabel(fClusterModuleLabel, clscol);
      art::fill_ptr_vector(clusters, clscol);
      this->CheckRecoClusters(clusters, scs, plist, allhits);
    }
    catch(cet::exception &e){
      mf::LogWarning("RecoCheckAna") << "could not recover clusters, failed with message:\n"
				     << e;
    }	  
  }
  if(fCheckTracks){
    try{
      e.getByLabel(fTrackModuleLabel, trkcol);
      art::fill_ptr_vector(tracks, trkcol);
      this->CheckRecoProngs(tracks, scs, plist, allhits, fTrackPurity, fTrackEfficiency);
    }
    catch(cet::exception &e){
      mf::LogWarning("RecoCheckAna") << "could not recover tracks, failed with message:\n"
				     << e;
    }	  
  }
  if(fCheckShowers){
    try{
      e.getByLabel(fShowerModuleLabel, shwcol);
      art::fill_ptr_vector(showers, shwcol);
      this->CheckRecoProngs(showers, scs, plist, allhits, fShowerPurity, fShowerEfficiency);
    }
    catch(cet::exception &e){
      mf::LogWarning("RecoCheckAna") << "could not recover showers, failed with message:\n"
				     << e;
    }	  
  }
  if(fCheckVertices){
    try{
      e.getByLabel(fVertexModuleLabel, vtxcol);
      art::fill_ptr_vector(vertices, vtxcol);
      this->CheckRecoVertices(vertices, scs, plist, allhits);
    }
    catch(cet::exception &e){
      mf::LogWarning("RecoCheckAna") << "could not recover vertices, failed with message:\n"
				     << e;
    }	  
  }
  if(fCheckEvents){
    try{
      e.getByLabel(fEventModuleLabel, evtcol);
      art::fill_ptr_vector(events, evtcol);
      this->CheckRecoEvents(events, scs, plist, allhits);
    }
    catch(cet::exception &e){
      mf::LogWarning("RecoCheckAna") << "could not recover events, failed with message:\n"
				     << e;
    }	  
  }

  return;
 
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::reconfigure(fhicl::ParameterSet const & p) 
{
  fHitModuleLabel     = p.get< std::string >("HitModuleLabel");
  fClusterModuleLabel = p.get< std::string >("ClusterModuleLabel");
  fShowerModuleLabel  = p.get< std::string >("ShowerModuleLabel"  );
  fTrackModuleLabel   = p.get< std::string >("TrackModuleLabel" );
  fVertexModuleLabel  = p.get< std::string >("VertexModuleLabel" );
  fEventModuleLabel   = p.get< std::string >("EventModuleLabel"  );
  fG4ModuleLabel      = p.get< std::string >("G4ModuleLabel"     );

  fCheckClusters      = p.get< bool        >("CheckClusters");
  fCheckShowers       = p.get< bool        >("CheckShowers" );
  fCheckTracks        = p.get< bool        >("CheckTracks"  );
  fCheckVertices      = p.get< bool        >("CheckVertices");
  fCheckEvents        = p.get< bool        >("CheckEvents"  );
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::beginRun(art::Run const &r) 
{
  art::ServiceHandle<art::TFileService> tfs;

  if(fCheckEvents){
    fEventPurity       = tfs->make<TH1D>("eventPurity",       ";Purity;Events",       100, 0., 1.1);
    fEventEfficiency   = tfs->make<TH1D>("eventEfficiency",   ";Efficiency;Events",   100, 0., 1.1);
  }
  if(fCheckVertices){
    fVertexPurity      = tfs->make<TH1D>("vertexPurity",      ";Purity;Vertices",     100, 0., 1.1);
    fVertexEfficiency  = tfs->make<TH1D>("vertexEfficiency",  ";Efficiency;Vertices", 100, 0., 1.1);
  }
  if(fCheckTracks){
    fTrackPurity       = tfs->make<TH1D>("trackPurity",       ";Purity;Tracks",       100, 0., 1.1);
    fTrackEfficiency   = tfs->make<TH1D>("trackEfficiency",   ";Efficiency;Tracks",   100, 0., 1.1);
  }
  if(fCheckShowers){
    fShowerPurity      = tfs->make<TH1D>("showerPurity",      ";Purity;Showers",      100, 0., 1.1);
    fShowerEfficiency  = tfs->make<TH1D>("showerEfficiency",  ";Efficiency;Showers",  100, 0., 1.1);
  }
  if(fCheckClusters){
    fClusterPurity     = tfs->make<TH1D>("clusterPurity",     ";Purity;Clusters",     110, 0., 1.1);
    fClusterEfficiency = tfs->make<TH1D>("clusterEfficiency", ";Efficiency;Clusters", 110, 0., 1.1);
  }
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoClusters(std::vector< art::Ptr<recob::Cluster> > const& clscol,
					    std::vector<const sim::SimChannel*>&           scs,
					    sim::ParticleList                       const& plist,
					    art::PtrVector<recob::Hit>              const& allhits)
{

  // grab the set of eve IDs for this plist
  std::set<int> eveIDs = cheat::BackTracker::GetSetOfEveIDs(plist);

  for(size_t c = 0; c < clscol.size(); ++c){

    // get the hits associated with this event
    art::PtrVector< recob::Hit > hits = clscol[c]->Hits();

    geo::View_t view = hits[0]->View();

    double maxPurity = -1.;
    double maxEfficiency = -1.;

    // loop over all eveIDs in this plist and determine which
    // has the highest purity and efficiency and put those values into the histogram
    
    std::set<int>::iterator itr = eveIDs.begin();
    while( itr != eveIDs.end() ){
    
      std::set<int> id;
      id.insert(*itr);
      // use the cheat::BackTracker to find purity and efficiency for these hits
      double purity     = cheat::BackTracker::HitCollectionPurity(plist, id, scs, hits);
      double efficiency = cheat::BackTracker::HitCollectionEfficiency(plist, id, scs, hits, allhits, view);

      // \todo Should the purity and efficiency be linked?  If so, then should the 
      // max purity be the deciding factor?
      if(purity     > maxPurity    ) maxPurity     = purity;
      if(efficiency > maxEfficiency) maxEfficiency = efficiency;

      itr++;
    }

    fClusterPurity    ->Fill(maxPurity);
    fClusterEfficiency->Fill(maxEfficiency);

  }// end loop over clusters

  return;
}

//-------------------------------------------------------------------
void cheat::RecoCheckAna::CheckRecoProngs(std::vector< art::Ptr<recob::Prong> > const& prgcol,
					  std::vector<const sim::SimChannel*>        & scs,
					  sim::ParticleList                     const& plist,
					  art::PtrVector<recob::Hit>            const& allhits,
					  TH1D*                                        purity,
					  TH1D*                                        efficiency)
{

  // grab the set of eve IDs for this plist
  std::set<int> eveIDs = cheat::BackTracker::GetSetOfEveIDs(plist);

  for(size_t p = 0; p < prgcol.size(); ++p){

    // get the hits associated with this event
    art::PtrVector< recob::Hit > hits = prgcol[p]->Hits();

    double maxPurity = -1.;
    double maxEfficiency = -1.;

    // loop over all eveIDs in this plist and determine which
    // has the highest purity and efficiency and put those values into the histogram
    
    std::set<int>::iterator itr = eveIDs.begin();
    while( itr != eveIDs.end() ){
    
      std::set<int> id;
      id.insert(*itr);
      // use the cheat::BackTracker to find purity and efficiency for these hits
      double purity     = cheat::BackTracker::HitCollectionPurity(plist, id, scs, hits);
      double efficiency = cheat::BackTracker::HitCollectionEfficiency(plist, id, scs, hits, allhits, geo::k3D);

      // \todo Should the purity and efficiency be linked?  If so, then should the 
      // max purity be the deciding factor?
      if(purity     > maxPurity    ) maxPurity     = purity;
      if(efficiency > maxEfficiency) maxEfficiency = efficiency;

      itr++;
    }

    purity    ->Fill(maxPurity);
    efficiency->Fill(maxEfficiency);

  }// end loop over events

  return;
}

//-------------------------------------------------------------------
//a true vertex will either consist of primary particles originating from
//the interaction vertex, or a primary particle decaying to make daughters
void cheat::RecoCheckAna::CheckRecoVertices(std::vector< art::Ptr<recob::Vertex> > const& vtxcol,
					    std::vector<const sim::SimChannel*>         & scs,
					    sim::ParticleList                      const& plist,
					    art::PtrVector<recob::Hit>             const& allhits)
{
  std::vector< std::set<int> > ids(1);
  // loop over all primary particles and put their ids into the first set of the 
  // vector.  add another set for each primary particle that also has daughters
  // and put those daughters into the new set
  for(int p = 0; p < plist.NumberOfPrimaries(); ++p){
    ids[0].insert(plist.Primary(p)->TrackId());
    if(plist.Primary(p)->NumberDaughters() > 0){
      std::set<int> dv;
      for(int d = 0; d < plist.Primary(p)->NumberDaughters(); ++d)
	dv.insert(plist.Primary(p)->Daughter(d));
      ids.push_back(dv);
    }//end if this primary particle has daughters
  }// end loop over primaries

  for(size_t v = 0; v < vtxcol.size(); ++v){

    // get the hits associated with this event
    art::PtrVector< recob::Hit > hits = vtxcol[v]->Hits();

    double maxPurity     = -1.;
    double maxEfficiency = -1.;

    for(size_t tv = 0; tv < ids.size(); ++tv){
      // use the cheat::BackTracker to find purity and efficiency for these hits
      double purity     = cheat::BackTracker::HitCollectionPurity(plist, ids[tv], scs, hits);
      double efficiency = cheat::BackTracker::HitCollectionEfficiency(plist, ids[tv], scs, hits, allhits, geo::k3D);

      if(purity     > maxPurity    ) maxPurity     = purity;
      if(efficiency > maxEfficiency) maxEfficiency = efficiency;
    }

    fVertexPurity    ->Fill(maxPurity);
    fVertexEfficiency->Fill(maxEfficiency);

  }// end loop over vertices

  return;
}

//-------------------------------------------------------------------
// in this method one should loop over the primary particles from a given
// MCTruth collection
// \todo need to divy it up in the case where there is more than 1 true interaction in a spill
void cheat::RecoCheckAna::CheckRecoEvents(std::vector< art::Ptr<recob::Event> > const& evtcol,
					  std::vector<const sim::SimChannel*>        & scs,
					  sim::ParticleList                     const& plist,
					  art::PtrVector<recob::Hit>            const& allhits)
{
  // loop over all primaries in the plist and grab them and their daughters to put into 
  // the set of track ids to pass on to the back tracker
  std::set<int> ids;
  for(int p = 0; p < plist.NumberOfPrimaries(); ++p){
    ids.insert(plist.Primary(p)->TrackId());
    for(int d = 0; d < plist.Primary(p)->NumberDaughters(); ++d)
      ids.insert(plist.Primary(p)->Daughter(d));
  }

  for(size_t ev = 0; ev < evtcol.size(); ++ev){

    // get the hits associated with this event
    art::PtrVector< recob::Hit > hits = evtcol[ev]->Hits();

    // use the cheat::BackTracker to find purity and efficiency for these hits
    double purity     = cheat::BackTracker::HitCollectionPurity(plist, ids, scs, hits);
    double efficiency = cheat::BackTracker::HitCollectionEfficiency(plist, ids, scs, hits, allhits, geo::k3D);

    fEventPurity    ->Fill(purity);
    fEventEfficiency->Fill(efficiency);

  }// end loop over events

  return;
}
