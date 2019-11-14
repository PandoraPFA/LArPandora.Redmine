/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Tracks
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace lar_pandora
{

    const std::vector<art::Ptr<recob::Hit>> LArPandoraTrackUtils::GetHits(const art::Ptr<recob::Track> &track, art::Event const &evt, const std::string &label)
    {    

        std::vector<art::Ptr<recob::Hit>> theseHits;
        GetAssocProductVector(track,evt,label,label,theseHits);
        return theseHits;
    }

    const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraTrackUtils::GetSpacePoints(const art::Ptr<recob::Track> &track, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::SpacePoint>> theseSPs;
        GetAssocProductVector(track,evt,label,label,theseSPs);
        return theseSPs;
    }

    const art::Ptr<recob::PFParticle> LArPandoraTrackUtils::GetParticle(const art::Ptr<recob::Track> &track, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::PFParticle>> theseParticles;
        GetAssocProductVector(track,evt,label,label,theseParticles);
        if (theseParticles.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraTrackUtils::GetParticle --- No associated particle found";
        }
        return theseParticles.at(0);
    }    

} // namespace lar_pandora


