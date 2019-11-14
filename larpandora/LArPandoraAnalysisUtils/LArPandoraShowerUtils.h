/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraShowerUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Showers
 *
 * @author leigh.howard.whitehead@cern.ch
*/

#ifndef LAR_PANDORA_SHOWER_UTILS_H
#define LAR_PANDORA_SHOWER_UTILS_H

#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

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

    static const std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Shower> &shower, art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Ptr<recob::Shower> &shower, art::Event const &evt, const std::string &label);

    static const art::Ptr<recob::PFParticle> GetParticle(const art::Ptr<recob::Shower> &shower, art::Event const &evt, const std::string &label);

private:

};

} // namespace lar_pandora


#endif // LAR_PANDORA_SHOWER_UTILS_H

