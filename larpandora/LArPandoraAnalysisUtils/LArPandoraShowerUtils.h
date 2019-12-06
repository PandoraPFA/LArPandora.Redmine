/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraShowerUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Showers
*/

#ifndef LAR_PANDORA_SHOWER_UTILS_H
#define LAR_PANDORA_SHOWER_UTILS_H

#include "art/Framework/Principal/Event.h"

// Access the type defs defined in the helper
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h"

#include <string>
#include <vector>

namespace lar_pandora
{

/**
 *
 * @brief LArPandoraShowerUtils class
 *
*/
class LArPandoraShowerUtils:LArPandoraUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the shower.
    *
    * @param shower is the shower for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the shower producer
    * 
    * @return vector of art::Ptrs to the hits
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the spacepoints associated with the shower.
    *
    * @param shower is the shower for which we want the spacepoints
    * @param evt is the underlying art event
    * @param label is the label for the shower producer
    * 
    * @return vector of art::Ptrs to the spacepoints
    */
    static std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the particle associated with the shower.
    *
    * @param shower is the shower for which we want the particle
    * @param evt is the underlying art event
    * @param label is the label for the shower producer
    * 
    * @return art::Ptr to the particle 
    */
    static art::Ptr<recob::PFParticle> GetPFParticle(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label);
};

} // namespace lar_pandora


#endif // LAR_PANDORA_SHOWER_UTILS_H

