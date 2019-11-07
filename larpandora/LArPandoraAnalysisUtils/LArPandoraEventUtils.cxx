/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.cxx
*
* @brief Utility containing helpful functions for end users to access products from events
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"

namespace lar_pandora
{

    const std::vector<art::Ptr<recob::PFParticle>> LArPandoraEventUtils::GetPFParticles(art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::PFParticle>> theseParticles;
        GetProductVector(evt,label,theseParticles);

        return theseParticles;
    }

    const std::vector<art::Ptr<recob::Track>> LArPandoraEventUtils::GetTracks(art::Event const &evt, const std::string &label)
    {
        mf::LogWarning("LArPandora") << " Please note: accessing PFParticle tracks through this method is not the recommended workflow.\n"
                                     << " Please use LArPandoraEventUtils::GetPFParticles and access the tracks with LArPandoraPFParticleUtils::GetTrack."
                                     << std::endl;

        std::vector<art::Ptr<recob::Track>> theseTracks;
        GetProductVector(evt,label,theseTracks);

        return theseTracks;
    }

    const std::vector<art::Ptr<recob::Shower>> LArPandoraEventUtils::GetShowers(art::Event const &evt, const std::string &label)
    {
        mf::LogWarning("LArPandora") << " Please note: accessing PFParticle showers through this method is not the recommended workflow.\n"
                                     << " Please use LArPandoraEventUtils::GetPFParticles and access the tracks with LArPandoraPFParticleUtils::GetShower."
                                     << std::endl;

        std::vector<art::Ptr<recob::Shower>> theseShowers;
        GetProductVector(evt,label,theseShowers);
     
        return theseShowers;
    }
} // namespace lar_pandora
