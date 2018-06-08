/**
 *  @file   larpandora/LArPandoraObjects/PFParticleMetadata.h
 *
 *  @brief  Metadata associated with a pandora produced PFParticle
 *
 */
#ifndef LAR_PANDORA_PFPARTICLE_METADATA_H
#define LAR_PANDORA_PFPARTICLE_METADATA_H

namespace larpandoraobj
{

/**
 *  @brief  Pandora PFParticleMetadata
 */
class PFParticleMetadata
{
public:
    /**
     *  @brief  Reconstruction hypothesis
     */
    enum ReconstructionHypothesis
    {
        COSMIC,
        NEUTRINO,
        TEST_BEAM
    };

    /**
     *  @brief  Constructor for clear cosmic rays
     */
    PFParticleMetadata();
    
    /**
     *  @brief  Constructor for particles in a slice
     *
     *  @param  hypothesis the reconstruction hypothesis used
     *  @param  sliceIndex the index of the slice in which the PFPartilce belongs
     */
    PFParticleMetadata(const ReconstructionHypothesis &hypothesis, const unsigned int sliceIndex);

    /**
     *  @brief  Get the reconstruction hypothesis
     *
     *  @return hypothesis
     */
    ReconstructionHypothesis GetReconstructionHypothesis() const;
    
    /**
     *  @brief  Check if the PFParticle is a clear cosmic ray
     *
     *  @return boolean
     */
    bool IsClearCosmicRay() const;

    /**
     *  @brief  Get the slice index
     *
     *  @return the slice index
     */
    unsigned int GetSliceIndex() const;

private:
    ReconstructionHypothesis m_hypothesis;    ///< The reconstruction hypothesis used
    bool                     m_isClearCosmic; ///< If the PFParticle is identified as a clear cosmic ray
    unsigned int             m_sliceIndex;    ///< The index of the slice (if available) in which the PFParticle belongs
};

} // namespace lar_pandora

#endif
