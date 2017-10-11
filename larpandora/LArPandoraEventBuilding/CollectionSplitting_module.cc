////////////////////////////////////////////////////////////////////////
// Class:       CollectionSplitting
// Plugin Type: producer (art v2_07_03)
// File:        CollectionSplitting_module.cc
//
// Generated at Wed Sep 27 12:37:53 2017 by Andrew D. Smith using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraEvent.h"

#include <memory>

namespace lar_pandora
{

class CollectionSplitting;


class CollectionSplitting : public art::EDProducer {
public:
  explicit CollectionSplitting(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CollectionSplitting(CollectionSplitting const &) = delete;
  CollectionSplitting(CollectionSplitting &&) = delete;
  CollectionSplitting & operator = (CollectionSplitting const &) = delete;
  CollectionSplitting & operator = (CollectionSplitting &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // FHicL congifurable parameters
  std::string     fInputProducerLabel;           ///< Label for the Pandora instance that produced the collections we want to split up
  std::string     fTrackProducerLabel;           ///< Label for the track producer using the Pandora instance that produced the collections we want to split up
  std::string     fShowerProducerLabel;          ///< Label for the shower producer using the Pandora instance that produced the collections we want to split up
  std::string     fHitProducerLabel;             ///< Label for the hit producer that was used as input to the Pandora instance specified
  bool            fShouldProduceNeutrinos;       ///< If we should produce collections related to neutrino top-level PFParticles
  bool            fShouldProduceCosmics;         ///< If we should produce collections related to cosmic top-level PFParticles
  bool            fShouldProduceT0s;             ///< If we should produce T0s (relevant when stitching over multiple drift volumes)
};


CollectionSplitting::CollectionSplitting(fhicl::ParameterSet const & p)
{

  reconfigure(p);

  // Define which types of collections this the module produces
  produces< std::vector<recob::PFParticle> >();
  produces< std::vector<recob::SpacePoint> >();
  produces< std::vector<recob::Cluster> >();

  if ( fShouldProduceT0s )
      produces< std::vector<anab::T0> >();

  produces< std::vector<recob::Vertex> >();
  produces< std::vector<recob::Track> >(); 
  produces< std::vector<recob::Shower> >();
  produces< std::vector<recob::PCAxis> >();

  produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
  produces< art::Assns<recob::PFParticle, recob::Cluster> >();

  if ( fShouldProduceT0s )
      produces< art::Assns<recob::PFParticle, anab::T0> >();

  produces< art::Assns<recob::PFParticle, recob::Vertex> >();
  produces< art::Assns<recob::PFParticle, recob::Track> >();
  produces< art::Assns<recob::PFParticle, recob::Shower> >();
  produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
  produces< art::Assns<recob::Track, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::PCAxis> >();
  produces< art::Assns<recob::SpacePoint, recob::Hit> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
}

void CollectionSplitting::produce(art::Event & e)
{
  if ( !fShouldProduceNeutrinos && !fShouldProduceCosmics ) 
    throw cet::exception("LArPandora") << " CollectionSplitting -- Must be configured to produce neutrinos or cosmics or both.";

  lar_pandora::LArPandoraEvent::Labels labels( fInputProducerLabel, fTrackProducerLabel, fShowerProducerLabel, fHitProducerLabel ); 
  lar_pandora::LArPandoraEvent fullEvent( this, &e, labels, fShouldProduceT0s );

  if ( fShouldProduceNeutrinos && fShouldProduceCosmics ) {
    fullEvent.WriteToEvent();
  }
  else {
    lar_pandora::LArPandoraEvent filteredEvent( fullEvent.FilterByPdgCode( fShouldProduceNeutrinos ) );
    filteredEvent.WriteToEvent();
  }
}

void CollectionSplitting::reconfigure(fhicl::ParameterSet const & p)
{
  fInputProducerLabel     = p.get<std::string>("InputProducerLabel");
  fTrackProducerLabel     = p.get<std::string>("TrackProducerLabel");
  fShowerProducerLabel    = p.get<std::string>("ShowerProducerLabel");
  fHitProducerLabel       = p.get<std::string>("HitProducerLabel");
  fShouldProduceNeutrinos = p.get<bool>("ShouldProduceNeutrinos", true);
  fShouldProduceCosmics   = p.get<bool>("ShouldProduceCosmics", true);
  fShouldProduceT0s       = p.get<bool>("ShouldProduceT0s", false);
}

DEFINE_ART_MODULE(CollectionSplitting)

} // namespace lar_pandora

