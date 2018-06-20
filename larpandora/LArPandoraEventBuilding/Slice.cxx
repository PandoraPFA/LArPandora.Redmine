/**
 *  @file   larpandora/LArPandoraEventBuilding/Slice.cxx
 *
 *  @brief  implementation of the lar pandora slice class
 */

#include "larpandora/LArPandoraEventBuilding/Slice.h"

namespace lar_pandora
{

Slice::Slice(const float nuScore, const PFParticleVector &nuHypothesis, const PFParticleVector &crHypothesis, const bool isNeutrino) : 
    m_nuScore(nuScore),
    m_nuHypothesis(nuHypothesis),
    m_crHypothesis(crHypothesis),
    m_isNeutrino(isNeutrino)
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

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool Slice::IsTaggedAsNeutrino() const
{
    return m_isNeutrino;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Slice::TagAsNeutrino()
{
    m_isNeutrino = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void Slice::TagAsCosmic()
{
    m_isNeutrino = false;
}

} // namespace lar_pandora
