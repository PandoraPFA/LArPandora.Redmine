/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraSlices.h
 *
 *  @brief  A mapping between two different reconstruction hypotheses on the same input hit collection
 */

#ifndef LAR_PANDORA_SLICES_H
#define LAR_PANDORA_SLICES_H 1

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Hit.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <memory>
#include <algorithm>
#include <map>

namespace lar_pandora
{

/**
 *  @brief LArPandoraSlices class
 */
class LArPandoraSlices
{
public:

    /**
     *  @brief  Constructor from an art::Event
     *
     *  @param  pProducer pointer to the producer to write the output
     *  @param  pEvent pointer to the event to process
     *  @param  crRecoProducerLabel label for the PFParticle producer using the CR hypothesis
     *  @param  nuRecoProducerLabel label for the PFParticle producer using the neutrino hypothesis
     *  @param  hitProducerLabel label for the hit producer used to make the input PFParticles
     */
    LArPandoraSlices(art::EDProducer *pProducer, art::Event *pEvent, const std::string &crRecoProducerLabel, const std::string &nuRecoProducerLabel,
        const std::string &hitProducerLabel);

    typedef unsigned int SliceId;

    /**
     *  @brief  Get the vector of slice Id
     */
    std::vector< SliceId > GetSlices();

    /**
     *  @brief  Get the vector of PFParticles in a given slice reconstructed as a cosmic ray
     *
     *  @param  sliceId the id of the slice to return
     */
    PFParticleVector GetSliceAsCR(const SliceId sliceId);

    /**
     *  @brief  Get the vector of PFParticles in a given slice reconstructed as a neutrino
     *
     *  @param  sliceId the id of the slice to return
     */
    PFParticleVector GetSliceAsNu(const SliceId sliceId);

    /**
     *  @brief  Identify a given slice as being the neutrino event
     *
     *  @param  sliceId the id of the slice to identify as containing the neutrino
     */
    void IdSliceAsNu(const SliceId sliceId);

    /**
     *  @brief  Produce the CR tags for the slices, and write to the event
     */
    void WriteTags();

private:
    /**
     *  @brief  Use the associations between PFParticles and hits to identify the slices
     */
    void IdentifySlices();

    /**
     *  @brief  Tag the PFParticles in the supplied vector as specified 
     *
     *  @param  shouldTagAsNeutrino should the particles be tagged as a neutrino (or cosmic)
     *  @param  pfParticleVector the input vector of PFParticles to be tagged 
     *  @param  outputTags the output collection of cosmic tags
     *  @param  outputAssn the output association from PFParticle -> CosmicTag
     */
    void WriteTag(const bool shouldTagAsNeutrino, const PFParticleVector &pfParticleVector, std::unique_ptr< std::vector< anab::CosmicTag > > &outputTags,
        std::unique_ptr< art::Assns< recob::PFParticle, anab::CosmicTag > > &outputAssn);

    /**
     *  @brief  Get a mapping between PFParticles and their Ids 
     *
     *  @param  inputParticles input PFParticles to obtain the mapping from
     *  @param  outputMap output mapping between PFParticles and their Ids
     */
    void GetPFParticleIdMap(const PFParticleVector &inputParticles, PFParticleMap &outputMap);
    
    /**
     *  @brief  Make a new slice for each top-level neutrino PFParticle supplied by filling m_nuSlicePFParticles. 
     *
     *  @param  nuPFParticleIdMap input mapping between PFParticles and their Ids
     *  @param  nuTopLevelParticles input vector of top-level PFParticles (i.e. the neutrino PFParticles)
     *  @param  nuFinalStateParticles input vector of final state PFParticles (i.e the first daughters of neutrino PFParticles)
     *  @param  nuFinalStateParticlesToSlice output mapping from final state PFParticles to slice id 
     */
    void MakeSlicePerNeutrino(const PFParticleMap &nuPFParticleIdMap, const PFParticleVector &nuTopLevelParticles, const PFParticleVector &nuFinalStateParticles,
        std::map< art::Ptr< recob::PFParticle>, SliceId > &nuFinalStateParticlesToSlice);

    /**
     *  @brief  Collect all PFParticles downstream of a supplied particle
     *
     *  @param  pfParticleMap input mapping between PFParticles and their Ids
     *  @param  part the particle to seed the collection
     *  @param  daughterParticles the output vector of all PFParticles downstream part
     */
    void CollectDaughters(const PFParticleMap &pfParticleMap, const art::Ptr< recob::PFParticle > &part, PFParticleVector &daughterParticles);

    /**
     *  @brief  Add CR PFParticles to existing neutrino slices (if they share hits), or make new slices
     *
     *  @param  crPFParticleIdMap  input mapping from cosmic ray PFParticles to their Ids
     *  @param  crFinalStateParticles input vector of final-state cosmic PFParticles 
     *  @param  crParticlesToHits input mapping from cosmic PFParticles to associated hits
     *  @param  nuHitsToParticles input mapping from hits to associated neutrino PFParticles
     *  @param  nuFinalStateParticlesToSlice input mapping from neutrino final-state PFParticles to slice ids
     */
    void AddCRParticlesToSlices(const PFParticleMap &crPFParticleIdMap, const PFParticleVector &crFinalStateParticles, const PFParticlesToHits &crParticlesToHits,
        const HitsToPFParticles &nuHitsToParticles, const std::map< art::Ptr< recob::PFParticle>, SliceId > &nuFinalStateParticlesToSlice);

    art::EDProducer    *m_pProducer;            ///< The producer which should write the output collections and associations
    art::Event         *m_pEvent;               ///< The event to consider

    const std::string   m_crRecoProducerLabel;  ///< Label for the PFParticle producer using the CR hypothesis
    const std::string   m_nuRecoProducerLabel;  ///< Label for the PFParticle producer using the neutrino hypothesis
    const std::string   m_hitProducerLabel;     ///< Label for the hit producer used to make the input PFParticles
    
    std::map< SliceId, PFParticleVector > m_crSlicePFParticles;  ///< Mapping from slice ID --> PFParticles reconstructed under the cosmic hypothesis
    std::map< SliceId, PFParticleVector > m_nuSlicePFParticles;  ///< Mapping from slice ID --> PFParticles reconstructed under the neutrino hypothesis

    bool                m_doesEventContainNeutrino;         ///< Has the user identified a slice as containing a neutrino
    SliceId             m_nuSliceId;                        ///< The identified slice that contains the neutrino (may not be set)
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_SLICES_H
