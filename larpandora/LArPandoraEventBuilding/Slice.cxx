/**
 *  @file   larpandora/LArPandoraEventBuilding/Slice.cxx
 *
 *  @brief  implementation of the lar pandora slice class
 */

#include "larpandora/LArPandoraEventBuilding/Slice.h"

namespace lar_pandora
{

Slice::Slice(const float nuScore, const PFParticleVector &nuHypothesis, const PFParticleVector &crHypothesis) : 
    m_nuScore(nuScore),
    m_nuHypothesis(nuHypothesis),
    m_crHypothesis(crHypothesis)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float Slice::GetNeutrinoScore() const
{
    return m_nuScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector Slice::GetNeutrinoHypothesis() const
{
    return m_nuHypothesis;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleVector Slice::GetCosmicRayHypothesis() const
{
    return m_crHypothesis;
}

} // namespace lar_pandora
