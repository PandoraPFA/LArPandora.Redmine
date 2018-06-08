/**
 *  @file   larpandora/LArPandoraObjects/PFParticleMetadata.cxx
 *
 *  @brief  Metadata associated with a pandora produced PFParticle implementation
 *
 */

#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

#include "cetlib_except/exception.h"

#include <limits>
#include <iostream>

namespace larpandoraobj
{

PFParticleMetadata::PFParticleMetadata() :
    m_hypothesis(COSMIC),
    m_isClearCosmic(true),
    m_sliceIndex(std::numeric_limits<unsigned int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleMetadata::PFParticleMetadata(const PFParticleMetadata::ReconstructionHypothesis &hypothesis, const unsigned int sliceIndex) :
    m_hypothesis(hypothesis),
    m_isClearCosmic(false),
    m_sliceIndex(sliceIndex)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
PFParticleMetadata::ReconstructionHypothesis PFParticleMetadata::GetReconstructionHypothesis() const
{
    return m_hypothesis;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

bool PFParticleMetadata::IsClearCosmicRay() const
{
    return m_isClearCosmic;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int PFParticleMetadata::GetSliceIndex() const
{
    if (this->IsClearCosmicRay())
        throw cet::exception("PFParticleMetadata") << " Can't access the slice index of a clear cosmic ray" << std::endl;

    return m_sliceIndex;
}

} // namespace lar_pandora
