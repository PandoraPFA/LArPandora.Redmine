////////////////////////////////////////////////////////////////////////
// Class:       LArPandoraNeutrinoId
// Plugin Type: producer (art v2_07_03)
// File:        LArPandoraNeutrinoId_module.cc
//
// Generated at Mon Oct  2 11:54:57 2017 by Andrew D. Smith using cetskelgen
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

#include "larpandora/LArPandoraEventBuilding/LArPandoraSlices.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraNeutrinoId;


class LArPandoraNeutrinoId : public art::EDProducer {
public:
  explicit LArPandoraNeutrinoId(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LArPandoraNeutrinoId(LArPandoraNeutrinoId const &) = delete;
  LArPandoraNeutrinoId(LArPandoraNeutrinoId &&) = delete;
  LArPandoraNeutrinoId & operator = (LArPandoraNeutrinoId const &) = delete;
  LArPandoraNeutrinoId & operator = (LArPandoraNeutrinoId &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.

  std::string   fCRProducerLabel;    ///<
  std::string   fNuProducerLabel;    ///<
  std::string   fHitProducerLabel;   ///<
};


LArPandoraNeutrinoId::LArPandoraNeutrinoId(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  reconfigure( p );

  // Call appropriate produces<>() functions here.
  produces< std::vector< anab::CosmicTag > >();
  produces< art::Assns< recob::PFParticle, anab::CosmicTag > >();
}

void LArPandoraNeutrinoId::produce(art::Event & e)
{
  // Implementation of required member function here.

  lar_pandora::LArPandoraSlices slices( this, &e, fCRProducerLabel, fNuProducerLabel, fHitProducerLabel );

  for ( const auto & slice : slices.GetSlices() ) {
    std::vector< art::Ptr< recob::PFParticle > > crParticles = slices.GetSliceAsCR( slice );
    std::vector< art::Ptr< recob::PFParticle > > nuParticles = slices.GetSliceAsNu( slice );

    // Logic to produce a map: sliceId → sliceProperties
    // ...
  }

  // Temporarily just choose the first slice.
  lar_pandora::LArPandoraSlices::SliceId bestSliceId = 0;

  slices.IdSliceAsNu( bestSliceId );
  slices.WriteTags();
}

void LArPandoraNeutrinoId::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fCRProducerLabel  = p.get<std::string>("CRProducerLabel");
  fNuProducerLabel  = p.get<std::string>("NuProducerLabel");
  fHitProducerLabel = p.get<std::string>("HitProducerLabel");
}

DEFINE_ART_MODULE(LArPandoraNeutrinoId)

} // namespace lar_pandora

