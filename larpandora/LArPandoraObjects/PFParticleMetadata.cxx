/**
 *  @file   larpandora/LArPandoraObjects/PFParticleMetadata.cxx
 *
 *  @brief  Metadata associated with a pandora produced PFParticle implementation
 *
 */

#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

namespace larpandoraobj
{

PFParticleMetadata::PFParticleMetadata()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleMetadata::PFParticleMetadata(const pandora::ParticleFlowObject *const pPfo) :
    m_propertiesMap(pPfo->GetPropertiesMap())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
const PFParticleMetadata::PropertiesMap &PFParticleMetadata::GetPropertiesMap() const
{
    return m_propertiesMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleMetadata::SetPropertiesMap(const PropertiesMap &propertiesMap)
{
    m_propertiesMap = propertiesMap;
}

} // namespace lar_pandora
