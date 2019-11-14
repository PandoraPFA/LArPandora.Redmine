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

    const std::vector<art::Ptr<recob::Hit>> LArPandoraShowerUtils::GetHits(const art::Ptr<recob::Shower> &shower, art::Event const &evt, const std::string &label)
    {    

        std::vector<art::Ptr<recob::Hit>> theseHits;
        GetAssocProductVector(shower,evt,label,label,theseHits);
        return theseHits;
    }

    const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraShowerUtils::GetSpacePoints(const art::Ptr<recob::Shower> &shower, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::SpacePoint>> theseSPs;
        GetAssocProductVector(shower,evt,label,label,theseSPs);
        return theseSPs;
    }

    const art::Ptr<recob::PFParticle> LArPandoraShowerUtils::GetParticle(const art::Ptr<recob::Shower> &shower, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::PFParticle>> theseParticles;
        GetAssocProductVector(shower,evt,label,label,theseParticles);
        if (theseParticles.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraShowerUtils::GetParticle --- No associated particle found";
        }
        return theseParticles.at(0);
    }    

} // namespace lar_pandora


