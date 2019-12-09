/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraSpacePointUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about SpacePoints
*/

#ifndef LAR_PANDORA_SPACEPOINT_UTILS_H
#define LAR_PANDORA_SPACEPOINT_UTILS_H

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
 * @brief LArPandoraSpacePointUtils class
 *
*/
class LArPandoraSpacePointUtils:LArPandoraUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the spacepoint.
    *
    * @param spacepoint is the spacepoint for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the spacepoint producer
    * 
    * @return vector of art::Ptrs to the hits 
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::SpacePoint> &pSpacepoint, const art::Event &evt, const std::string &label);
};

} // namespace lar_pandora


#endif // LAR_PANDORA_SPACEPOINT_UTILS_H

