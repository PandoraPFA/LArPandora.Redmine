/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
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

    const std::vector<art::Ptr<anab::T0>> LArPandoraPFParticleUtils::GetT0(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<anab::T0>> theseT0s;
        GetAssocProductVector(part,evt,label,label,theseT0s);
        return theseT0s;
    }

    const std::vector<art::Ptr<anab::CosmicTag>> LArPandoraPFParticleUtils::GetCosmicTag(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<anab::CosmicTag>> theseTags;
        GetAssocProductVector(part,evt,label,label,theseTags); 
        return theseTags;
    }

    const std::vector<art::Ptr<recob::PFParticle>> LArPandoraPFParticleUtils::GetDaughterParticles(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

        art::Handle<std::vector<recob::PFParticle>> particles;
        evt.getByLabel(label,particles);

        unsigned int thisParticleIndex = part.key();

        std::vector<art::Ptr<recob::PFParticle>> daughters;

        for (unsigned int p = 0; p < particles->size(); ++p)
        {     
            if (particles->at(p).Parent() == thisParticleIndex)
            {
                art::Ptr<recob::PFParticle> pPtr(particles,p);
            }
        }

        return daughters;
    }

    const std::vector<art::Ptr<recob::Hit>> LArPandoraPFParticleUtils::GetHits(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {    

        std::vector<art::Ptr<recob::Hit>> theseHits;
        GetAssocProductVector(part,evt,label,label,theseHits);
        return theseHits;
    }

    const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraPFParticleUtils::GetSpacePoints(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::SpacePoint>> theseSPs;
        GetAssocProductVector(part,evt,label,label,theseSPs);
        return theseSPs;
    }

    const art::Ptr<recob::Track> LArPandoraPFParticleUtils::GetTrack(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &particleLabel, const std::string &trackLabel)
    {

        std::vector<art::Ptr<recob::Track>> theseTracks;
        GetAssocProductVector(part,evt,particleLabel,trackLabel,theseTracks);
        if (theseTracks.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraPFParticleUtils::GetTrack --- No associated track found";
        }
        return theseTracks.at(0);
    }    
    
    const art::Ptr<recob::Shower> LArPandoraPFParticleUtils::GetShower(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &particleLabel, const std::string &showerLabel)
    {

        std::vector<art::Ptr<recob::Shower>> theseShowers;
        GetAssocProductVector(part,evt,particleLabel,showerLabel,theseShowers);

        if (theseShowers.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraPFParticleUtils::GetShower --- No associated shower found";
        }
        return theseShowers.at(0);

    }

    const art::Ptr<larpandoraobj::PFParticleMetadata> LArPandoraPFParticleUtils::GetMetadata(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> theseMetadata;
        GetAssocProductVector(part,evt,label,label,theseMetadata);

        if (theseMetadata.size() == 0)
        {
            throw cet::exception("LArPandora") << " LArPandoraPFParticleUtils::GetMetadata --- No associated metadata found";
        }
        return theseMetadata.at(0);
    }

    bool LArPandoraPFParticleUtils::IsTrack(const art::Ptr<recob::PFParticle> particle)
    {
        return LArPandoraHelper::IsTrack(particle);
    }

    bool LArPandoraPFParticleUtils::IsShower(const art::Ptr<recob::PFParticle> particle)
    {
        return LArPandoraHelper::IsShower(particle);
    }


} // namespace lar_pandora


