/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.cxx
*
* @brief Utility containing helpful functions for end users to access products from events
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::PFParticle>> LArPandoraEventUtils::GetPFParticles(const art::Event &evt, const std::string &label)
{
    return LArPandoraEventUtils::GetProductVector<recob::PFParticle>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Track>> LArPandoraEventUtils::GetTracks(const art::Event &evt, const std::string &label)
{
    mf::LogWarning("LArPandora") << " Please note: accessing PFParticle tracks through this method is not the recommended workflow.\n"
                                 << " Please use LArPandoraEventUtils::GetPFParticles and access the tracks with LArPandoraPFParticleUtils::GetTrack."
                                 << std::endl;

    return LArPandoraEventUtils::GetProductVector<recob::Track>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Shower>> LArPandoraEventUtils::GetShowers(const art::Event &evt, const std::string &label)
{
    mf::LogWarning("LArPandora") << " Please note: accessing PFParticle showers through this method is not the recommended workflow.\n"
                                 << " Please use LArPandoraEventUtils::GetPFParticles and access the tracks with LArPandoraPFParticleUtils::GetShower."
                                 << std::endl;

    return LArPandoraEventUtils::GetProductVector<recob::Shower>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Vertex>> LArPandoraEventUtils::GetVertices(const art::Event &evt, const std::string &label)
{
    return LArPandoraEventUtils::GetProductVector<recob::Vertex>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> LArPandoraEventUtils::GetSpacePoints(const art::Event &evt, const std::string &label)
{
    return LArPandoraEventUtils::GetProductVector<recob::SpacePoint>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Slice>> LArPandoraEventUtils::GetSlices(const art::Event &evt, const std::string &label)
{
    return LArPandoraEventUtils::GetProductVector<recob::Slice>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::PFParticle>> LArPandoraEventUtils::GetClearCosmics(const art::Event &evt, const std::string &label)
{
    std::vector<art::Ptr<recob::PFParticle>> theseParticles = LArPandoraEventUtils::GetProductVector<recob::PFParticle>(evt,label);

    std::vector<art::Ptr<recob::PFParticle>> theseCosmics;

    // We only want primary cosmic rays
    for (art::Ptr<recob::PFParticle> pParticle : theseParticles)
    {
        if (!pParticle->IsPrimary())
        {
            continue;
        } 

        if (LArPandoraPFParticleUtils::IsClearCosmic(pParticle, evt, label))
        {
            theseCosmics.push_back(pParticle);
        }
    }

    return theseCosmics;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> LArPandoraEventUtils::GetNeutrino(const art::Event &evt, const std::string &label)
{
    if (!HasNeutrino(evt,label))
    {
        throw cet::exception("LArPandora") << "LArPandoraEventUtils::GetNeutrino --- No neutrino found";
    }
    
    art::Ptr<recob::PFParticle> pNeutrino;
    std::vector<art::Ptr<recob::PFParticle>> particles = LArPandoraEventUtils::GetPFParticles(evt,label);
    for (art::Ptr<recob::PFParticle> pParticle : particles)
    {
        if (LArPandoraPFParticleUtils::IsNeutrino(pParticle))
        {
            pNeutrino = pParticle;
            break; 
        }
    }
    return pNeutrino;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraEventUtils::HasNeutrino(const art::Event &evt, const std::string &label)
{
    bool hasNeutrino = false;
    std::vector<art::Ptr<recob::PFParticle>> particles = LArPandoraEventUtils::GetPFParticles(evt,label);
    for (art::Ptr<recob::PFParticle> pParticle : particles)
    {
        if(!pParticle->IsPrimary())
        {
            continue;
        }
        if (LArPandoraPFParticleUtils::IsNeutrino(pParticle))
        {
            hasNeutrino = true;
            break;
        }
    }
    return hasNeutrino;
}

} // namespace lar_pandora

