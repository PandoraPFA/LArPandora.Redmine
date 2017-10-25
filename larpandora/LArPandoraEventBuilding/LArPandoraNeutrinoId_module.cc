/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraNeutrinoId.cc
 *
 *  @brief  module for lar pandora neutrino identification
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraNeutrinoId : public art::EDProducer {
public:
    explicit LArPandoraNeutrinoId(fhicl::ParameterSet const & pset);

    LArPandoraNeutrinoId(LArPandoraNeutrinoId const &) = delete;
    LArPandoraNeutrinoId(LArPandoraNeutrinoId &&) = delete;
    LArPandoraNeutrinoId & operator = (LArPandoraNeutrinoId const &) = delete;
    LArPandoraNeutrinoId & operator = (LArPandoraNeutrinoId &&) = delete;

    void produce(art::Event & evt) override;

private:
    std::string  m_CRProducerLabel;    ///<  the input PFParticle producer under the cosmic hypothesis
    std::string  m_NuProducerLabel;    ///<  the input PFParticle producer under the neutrino hypothesis
    std::string  m_HitProducerLabel;   ///<  the input hit producer used to produces the input PFParticles
};

DEFINE_ART_MODULE(LArPandoraNeutrinoId)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraSlices.h"

namespace lar_pandora
{

LArPandoraNeutrinoId::LArPandoraNeutrinoId(fhicl::ParameterSet const & pset) :
    m_CRProducerLabel(pset.get<std::string>("CRProducerLabel")),
    m_NuProducerLabel(pset.get<std::string>("NuProducerLabel")),
    m_HitProducerLabel(pset.get<std::string>("HitProducerLabel"))
{
    produces< std::vector< anab::CosmicTag > >();
    produces< art::Assns< recob::PFParticle, anab::CosmicTag > >();
}

void LArPandoraNeutrinoId::produce(art::Event & evt)
{
    lar_pandora::LArPandoraSlices slices(this, &evt, m_CRProducerLabel, m_NuProducerLabel, m_HitProducerLabel);

    for (const auto & slice : slices.GetSlices()) {
        std::vector< art::Ptr< recob::PFParticle > > crParticles = slices.GetSliceAsCR(slice);
        std::vector< art::Ptr< recob::PFParticle > > nuParticles = slices.GetSliceAsNu(slice);

        // Logic to produce a map: sliceId --> sliceProperties
        // ...
    }

    // Temporarily just choose the first slice.
    lar_pandora::LArPandoraSlices::SliceId bestSliceId = 0;

    slices.IdSliceAsNu(bestSliceId);
    slices.WriteTags();
}

} // namespace lar_pandora

