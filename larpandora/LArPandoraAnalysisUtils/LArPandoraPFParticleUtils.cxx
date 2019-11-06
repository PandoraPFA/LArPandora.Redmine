/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

namespace lar_pandora
{

    const std::vector<art::Ptr<anab::T0>> LArPandoraPFParticleUtils::GetT0(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<anab::T0>> theseT0s;
        GetAssocProductVector(part,evt,label,label,theseT0s);
//        art::Handle<std::vector<recob::PFParticle>> particles;
//        evt.getByLabel(label,particles);
//
//        if (!particles.isValid())
//        {
//            mf::LogError("LArPandora") << " Failed to find PFParticles... returning empty vector" << std::endl;
//            return std::vector<art::Ptr<anab::T0>>();
//        }
//
//        const art::FindManyP<anab::T0> findParticleT0s(particles,evt,label);
//    
//        const std::vector<art::Ptr<anab::T0>> theseT0s = findParticleT0s.at(part.key());
//
        return theseT0s;
    }

    const std::vector<art::Ptr<anab::CosmicTag>> LArPandoraPFParticleUtils::GetCosmicTag(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

          std::vector<art::Ptr<anab::CosmicTag>> theseTags;
          GetAssocProductVector(part,evt,label,label,theseTags); 

//        art::Handle<std::vector<recob::PFParticle>> particles;
//        evt.getByLabel(label,particles);
//
//        if (!particles.isValid())
//        {
//            mf::LogError("LArPandora") << " Failed to find PFParticles... returning empty vector" << std::endl;
//            return std::vector<art::Ptr<anab::CosmicTag>>();
//        }
//
//        const art::FindManyP<anab::CosmicTag> findParticleTags(particles,evt,label);
//    
//        const std::vector<art::Ptr<anab::CosmicTag>> theseTags = findParticleTags.at(part.key());

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
//        art::Handle<std::vector<recob::PFParticle>> particles;
//        evt.getByLabel(label,particles);
//        
//        if (!particles.isValid())
//        {   
//            mf::LogError("LArPandora") << " Failed to find PFParticles... returning empty vector" << std::endl;
//            return std::vector<art::Ptr<recob::Hit>>();
//        }
//        
//        const art::FindManyP<recob::Hit> findParticleHits(particles,evt,label);
//        
//        const std::vector<art::Ptr<recob::Hit>> theseHits = findParticleHits.at(part.key());
//        
        return theseHits;
    }

    const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraPFParticleUtils::GetSpacePoints(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::SpacePoint>> theseSPs;
        GetAssocProductVector(part,evt,label,label,theseSPs);
//        art::Handle<std::vector<recob::PFParticle>> particles;
//        evt.getByLabel(label,particles);
//        
//        if (!particles.isValid())
//        {   
//            mf::LogError("LArPandora") << " Failed to find PFParticles... returning empty vector" << std::endl;
//            return std::vector<art::Ptr<recob::SpacePoint>>();
//        }
//        
//        const art::FindManyP<recob::SpacePoint> findParticleSPs(particles,evt,label);
//        
//        const std::vector<art::Ptr<recob::SpacePoint>> theseSPs = findParticleSPs.at(part.key());
        
        return theseSPs;
    }

    const art::Ptr<recob::Track> LArPandoraPFParticleUtils::GetTrack(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &particleLabel, const std::string &trackLabel)
    {

        std::vector<art::Ptr<recob::Track>> theseTracks;
        GetAssocProductVector(part,evt,particleLabel,trackLabel,theseTracks);
//        // Pandora produces associations between PFParticles and recob::Track objects
//        auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(label);
//        const art::FindManyP<recob::Track> findTracks(particles,evt,label);
//        const std::vector<art::Ptr<recob::Track>> pfpTracks = findTracks.at(part->Self());
//        // Check that the track exists
//        if(pfpTracks.size() == 0){
//          mf::LogError("LArPandora") << "This particle has no track. Returning nullptr." << std::endl;
//          assert(0);
//        }
        if (theseTracks.size() == 0)
        {
          mf::LogError("LArPandora") << "No associated track found... exiting." << std::endl;
          assert(0);
        }
        return theseTracks.at(0);
    }    
    
    const art::Ptr<recob::Shower> LArPandoraPFParticleUtils::GetShower(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &particleLabel, const std::string &showerLabel)
    {

        std::vector<art::Ptr<recob::Shower>> theseShowers;
        GetAssocProductVector(part,evt,particleLabel,showerLabel,theseShowers);

//        // Pandora produces associations between PFParticles and recob::Track objects
//        auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(label);
//        const art::FindManyP<recob::Shower> findShowers(particles,evt,label);
//        const std::vector<art::Ptr<recob::Shower>> pfpShowers = findShowers.at(part->Self());
//        // Check that the track exists
//        if(pfpShowers.size() == 0){
//          mf::LogError("LArPandora") << "This particle has no track. Returning nullptr." << std::endl;
//          assert(0);
//        }

        if (theseShowers.size() == 0)
        {
          mf::LogError("LArPandora") << "No associated shower found... exiting." << std::endl;
          assert(0);
        }
        return theseShowers.at(0);

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


