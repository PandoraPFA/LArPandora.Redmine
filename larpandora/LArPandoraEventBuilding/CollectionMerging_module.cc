////////////////////////////////////////////////////////////////////////
// Class:       CollectionMerging
// Plugin Type: producer (art v2_07_03)
// File:        CollectionMerging_module.cc
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
    
class CollectionMerging;


class CollectionMerging : public art::EDProducer {
public:
  explicit CollectionMerging(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CollectionMerging(CollectionMerging const &) = delete;
  CollectionMerging(CollectionMerging &&) = delete;
  CollectionMerging & operator = (CollectionMerging const &) = delete;
  CollectionMerging & operator = (CollectionMerging &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // FHicL congifurable parameters
  std::string  fAllHitsCRProducerLabel;          ///< Label of the pandora instance that ran CR reco on all hits
  std::string  fAllHitsCRTrackProducerLabel;     ///< Label of the track producer using the pandora instance that ran CR reco on all hits
  std::string  fAllHitsCRShowerProducerLabel;    ///< Label of the shower producer using the pandora instance that ran CR reco on all hits

  std::string  fCRRemHitsCRProducerLabel;        ///< Label of the pandora instance that ran CR reco on CR removed hits
  std::string  fCRRemHitsCRTrackProducerLabel;   ///< Label of the track producer using the pandora instance that ran CR reco on CR removed hits
  std::string  fCRRemHitsCRShowerProducerLabel;  ///< Label of the shower producer using the pandora instance that ran CR reco on CR removed hits

  std::string  fCRRemHitsNuProducerLabel;        ///< Label of the pandora instance that ran Nu reco on CR removed hits
  std::string  fCRRemHitsNuTrackProducerLabel;   ///< Label of the track producer using the pandora instance that ran Nu reco on CR removed hits
  std::string  fCRRemHitsNuShowerProducerLabel;  ///< Label of the shower producer using the pandora instance that ran Nu reco on CR removed hits

  std::string  fAllHitProducerLabel;             ///< Label of the primary hit producer 
  std::string  fCRRemHitProducerLabel;           ///< Label of the CR removed hit producer

  std::string  fClearCRTagProducerLabel;         ///< Label of the unabiguous CR tag producer
  std::string  fNuIdCRTagProducerLabel;          ///< Label of the neutrino-ID CR tag producer

  bool         fShouldProduceNeutrinos;          ///< If we should produce collections related to neutrino top-level PFParticles
  bool         fShouldProduceT0s;                ///< If we should produce T0s (relevant when stitching over multiple drift volumes)
};


CollectionMerging::CollectionMerging(fhicl::ParameterSet const & p)
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

void CollectionMerging::produce(art::Event & e)
{
  // Get all reconstructions of the event
  lar_pandora::LArPandoraEvent::Labels allHitsCRLabels( fAllHitsCRProducerLabel, fAllHitsCRTrackProducerLabel, fAllHitsCRShowerProducerLabel, fAllHitProducerLabel );
  lar_pandora::LArPandoraEvent allHitsCREvent(   this, &e, allHitsCRLabels, fShouldProduceT0s   );

  lar_pandora::LArPandoraEvent::Labels crRemHitsCRLabels( fCRRemHitsCRProducerLabel, fCRRemHitsCRTrackProducerLabel, fCRRemHitsCRShowerProducerLabel, fCRRemHitProducerLabel );
  lar_pandora::LArPandoraEvent crRemHitsCREvent( this, &e, crRemHitsCRLabels, fShouldProduceT0s );

  lar_pandora::LArPandoraEvent::Labels crRemHitsNuLabels( fCRRemHitsNuProducerLabel, fCRRemHitsNuTrackProducerLabel, fCRRemHitsNuShowerProducerLabel, fCRRemHitProducerLabel );
  lar_pandora::LArPandoraEvent crRemHitsNuEvent( this, &e, crRemHitsNuLabels, fShouldProduceT0s );

  // Filter and merge into a consolidated output
  if ( fShouldProduceNeutrinos ) {
    lar_pandora::LArPandoraEvent filteredCRRemHitsNuEvent( crRemHitsNuEvent.FilterByCRTag( fShouldProduceNeutrinos, fNuIdCRTagProducerLabel  ) );
    filteredCRRemHitsNuEvent.WriteToEvent();
  }
  else{
    lar_pandora::LArPandoraEvent filteredAllHitsCREvent(   allHitsCREvent.FilterByCRTag(   fShouldProduceNeutrinos, fClearCRTagProducerLabel ) );
    lar_pandora::LArPandoraEvent filteredCRRemHitsCREvent( crRemHitsCREvent.FilterByCRTag( fShouldProduceNeutrinos, fNuIdCRTagProducerLabel  ) );
    lar_pandora::LArPandoraEvent mergedEvent( filteredAllHitsCREvent.Merge( filteredCRRemHitsCREvent ) );
    mergedEvent.WriteToEvent();
  }
}

void CollectionMerging::reconfigure(fhicl::ParameterSet const & p)
{
  fAllHitsCRProducerLabel         = p.get<std::string>("AllHitsCRProducerLabel");
  fAllHitsCRTrackProducerLabel    = p.get<std::string>("AllHitsCRTrackProducerLabel");
  fAllHitsCRShowerProducerLabel   = p.get<std::string>("AllHitsCRShowerProducerLabel");

  fCRRemHitsCRProducerLabel       = p.get<std::string>("CRRemHitsCRProducerLabel");
  fCRRemHitsCRTrackProducerLabel  = p.get<std::string>("CRRemHitsCRTrackProducerLabel");
  fCRRemHitsCRShowerProducerLabel = p.get<std::string>("CRRemHitsCRShowerProducerLabel");

  fCRRemHitsNuProducerLabel       = p.get<std::string>("CRRemHitsNuProducerLabel");
  fCRRemHitsNuTrackProducerLabel  = p.get<std::string>("CRRemHitsNuTrackProducerLabel");
  fCRRemHitsNuShowerProducerLabel = p.get<std::string>("CRRemHitsNuShowerProducerLabel");

  fAllHitProducerLabel      = p.get<std::string>("AllHitProducerLabel");
  fCRRemHitProducerLabel    = p.get<std::string>("CRRemHitProducerLabel");

  fClearCRTagProducerLabel  = p.get<std::string>("ClearCRTagProducerLabel");
  fNuIdCRTagProducerLabel   = p.get<std::string>("NuIdCRTagProducerLabel");

  fShouldProduceNeutrinos   = p.get<bool>("ShouldProduceNeutrinos", true);
  fShouldProduceT0s         = p.get<bool>("ShouldProduceT0s", false);
}

DEFINE_ART_MODULE(CollectionMerging)

} // namespace lar_pandora

