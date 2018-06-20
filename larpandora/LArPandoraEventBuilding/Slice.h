/**
 *  @file   larpandora/LArPandoraEventBuilding/Slice.h
 *
 *  @brief  header for the lar pandora slice class
 */

#ifndef LAR_PANDORA_SLICE_H
#define LAR_PANDORA_SLICE_H 1

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace lar_pandora
{

/**
 *  @brief Slice class
 */
class Slice
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  nuScore the neutrino score from Pandora
     *  @param  nuHypothesis the slice as reconstructed under the neutrino hypothesis
     *  @param  crHypothesis the slice as reconstructed under the cosmic-ray hypothesis
     */
    Slice(const float nuScore, const PFParticleVector &nuHypothesis, const PFParticleVector &crHypothesis);

    /**
     *  @brief Get the neutrino score for the slice
     */
    float GetNeutrinoScore() const;

    /**
     *  @brief Get the slice as reconstructed under the neutrino hypothesis
     */
    PFParticleVector GetNeutrinoHypothesis() const;
    
    /**
     *  @brief Get the slice as reconstructed under the cosmic-ray hypothesis
     */
    PFParticleVector GetCosmicRayHypothesis() const;

private:
    float m_nuScore;                    ///< The neutrino score from Pandora
    PFParticleVector m_nuHypothesis;    ///< The slice as reconstructed under the neutrino hypothesis
    PFParticleVector m_crHypothesis;    ///< The slice as reconstructed under the cosmic-ray hypothesis
};

typedef std::vector<Slice> SliceVector;

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_SLICE_H
