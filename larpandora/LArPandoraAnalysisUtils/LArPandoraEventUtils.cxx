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

    const std::vector<art::Ptr<recob::Vertex>> LArPandoraEventUtils::GetVertices(art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::Vertex>> theseVertices;
        GetProductVector(evt,label,theseVertices);

        return theseVertices;
    }

    const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraEventUtils::GetSpacePoints(art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::SpacePoint>> theseSpacePoints;
        GetProductVector(evt,label,theseSpacePoints);

        return theseSpacePoints;
    }

    const std::vector<art::Ptr<recob::Slice>> LArPandoraEventUtils::GetSlices(art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::Slice>> theseSlices;
        GetProductVector(evt,label,theseSlices);

        return theseSlices;
    }

    const std::vector<art::Ptr<recob::PFParticle>> LArPandoraEventUtils::GetClearCosmics(art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::PFParticle>> theseParticles;
        GetProductVector(evt,label,theseParticles);

        std::vector<art::Ptr<recob::PFParticle>> theseCosmics;

        for (art::Ptr<recob::PFParticle> particle : theseParticles)
        {
            if(LArPandoraPFParticleUtils::IsClearCosmic(particle, evt, label))
            {
                theseCosmics.push_back(particle);
            }
        }

        return theseCosmics;
    }

    const art::Ptr<recob::PFParticle> LArPandoraEventUtils::GetNeutrino(art::Event const &evt, const std::string &label)
    {

        if (!HasNeutrino(evt,label))
        {
          throw cet::exception("LArPandora") << "LArPandoraEventUtils::GetNeutrino --- No neutrino found";
        }
        
        art::Ptr<recob::PFParticle> neutrino;
        const std::vector<art::Ptr<recob::PFParticle>> particles = GetPFParticles(evt,label);
        for (art::Ptr<recob::PFParticle> particle : particles)
        {
            if (LArPandoraPFParticleUtils::IsNeutrino(particle))
            {
                neutrino = particle;      
            }
        }
        return neutrino;
    }

    const bool LArPandoraEventUtils::HasNeutrino(art::Event const &evt, const std::string &label)
    {
        bool hasNeutrino = false;
        const std::vector<art::Ptr<recob::PFParticle>> particles = GetPFParticles(evt,label);
        for (art::Ptr<recob::PFParticle> particle : particles)
        {
            if (LArPandoraPFParticleUtils::IsNeutrino(particle))
            {
                hasNeutrino = true;
                break;
            }
        }
        return hasNeutrino;
    }

} // namespace lar_pandora

