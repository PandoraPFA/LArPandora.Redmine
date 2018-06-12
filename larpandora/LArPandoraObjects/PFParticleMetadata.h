/**
 *  @file   larpandora/LArPandoraObjects/PFParticleMetadata.h
 *
 *  @brief  Metadata associated with a pandora produced PFParticle
 *
 */
#ifndef LAR_PANDORA_PFPARTICLE_METADATA_H
#define LAR_PANDORA_PFPARTICLE_METADATA_H

#include "Objects/ParticleFlowObject.h"

#include <map>
#include <string>

namespace larpandoraobj
{

/**
 *  @brief  Pandora PFParticleMetadata
 */
class PFParticleMetadata
{
public:
    /**
     *  @brief  Default constructor
     */
    PFParticleMetadata();

    /**
     *  @brief  Constructor
     *
     *  @param  hypothesis the reconstruction hypothesis used
     */
    PFParticleMetadata(const pandora::ParticleFlowObject *const pPfo);

    typedef pandora::PropertiesMap PropertiesMap;   ///< The properties map typedef

    /**
     *  @brief  Get the properties map
     *
     *  @return const reference to the properties map
     */
    const PropertiesMap &GetPropertiesMap() const;

    /**
     *  @brief  Set the properties map, or replace the existing map, with a user-provided map
     *
     *  @param  propertiesMap the replacement properties map
     */
    void SetPropertiesMap(const PropertiesMap &propertiesMap);

    // Further structures/properties to be added as they are requested

private:
    PropertiesMap   m_propertiesMap;                ///< The properties map
};

} // namespace lar_pandora

#endif
