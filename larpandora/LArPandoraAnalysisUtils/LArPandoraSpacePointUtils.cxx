/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraSpacePointUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about SpacePoints
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraSpacePointUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::Hit>> LArPandoraSpacePointUtils::GetHits(const art::Ptr<recob::SpacePoint> &pSpacepoint, const art::Event &evt, const std::string &label)
{    
    return LArPandoraSpacePointUtils::GetAssocProductVector<recob::Hit>(pSpacepoint,evt,label,label);
}

} // namespace lar_pandora

