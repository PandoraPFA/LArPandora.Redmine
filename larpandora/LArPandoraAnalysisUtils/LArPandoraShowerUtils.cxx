/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraShowerUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Showers
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraShowerUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::Hit>> LArPandoraShowerUtils::GetHits(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label)
{    
    return LArPandoraShowerUtils::GetAssocProductVector<recob::Hit>(pShower,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> LArPandoraShowerUtils::GetSpacePoints(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label)
{
    return LArPandoraShowerUtils::GetAssocProductVector<recob::SpacePoint>(pShower,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> LArPandoraShowerUtils::GetPFParticle(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label)
{
    return LArPandoraShowerUtils::GetAssocProduct<recob::PFParticle>(pShower,evt,label,label);
}    

} // namespace lar_pandora


