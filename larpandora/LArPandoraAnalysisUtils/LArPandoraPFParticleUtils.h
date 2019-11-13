/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about PFParticles
 *
 * @author leigh.howard.whitehead@cern.ch
*/

#ifndef LAR_PANDORA_PFPARTICLE_UTILS_H
#define LAR_PANDORA_PFPARTICLE_UTILS_H

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
 * @brief LArPandoraPFParticleUtils class
 *
*/
class LArPandoraPFParticleUtils:LArPandoraUtilsBase
{

public:

    /**
    * @brief Get the T0(s) associated with the particle.
    *
    * @param part particle for which we want the T0
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    */ 
    static const std::vector<art::Ptr<anab::T0>> GetT0(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label);

    /**
    * @brief Get the Cosmic Tag(s) associated with the particle.
    *
    * @param part particle for which we want the cosmic tag
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    */
    static const std::vector<art::Ptr<anab::CosmicTag>> GetCosmicTag(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::PFParticle>> GetDaughterParticles(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label);

    static const std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label);

    static const art::Ptr<recob::Track> GetTrack(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &particleLabel, const std::string &trackLabel);

    static const art::Ptr<recob::Shower> GetShower(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &particleLabel, const std::string &showerLabel);

    static const art::Ptr<larpandoraobj::PFParticleMetadata> GetMetadata(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label);

    static bool IsTrack(const art::Ptr<recob::PFParticle> particle);

    static bool IsShower(const art::Ptr<recob::PFParticle> particle);

    static bool IsNeutrino(const art::Ptr<recob::PFParticle particle);

private:

//    template <typename T> static void GetAssocProductVector(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label, const std::string &assocLabel, std::vector<art::Ptr<T>> &productVector);

};

    // Implementation of the template function to get the associated products from the event
//    template <typename T> void LArPandoraPFParticleUtils::GetAssocProductVector(const art::Ptr<recob::PFParticle> part, art::Event const &evt, const std::string &label, const std::string &assocLabel, std::vector<art::Ptr<T>> &productVector)
//    {
//
//        art::Handle<std::vector<recob::PFParticle>> particles;
//        evt.getByLabel(label,particles);
//
//        if (!particles.isValid())
//        {
//            mf::LogError("LArPandora") << " Failed to find PFParticles... returning empty vector" << std::endl;
//            productVector = std::vector<art::Ptr<T>>();
//            return;
//        }
//
//        const art::FindManyP<T> findParticleAssocs(particles,evt,assocLabel);
//
//        const std::vector<art::Ptr<T>> theseAssocs = findParticleAssocs.at(part.key());
//
//        productVector = theseAssocs;
//    }


} // namespace lar_pandora


#endif // LAR_PANDORA_PFPARTICLE_UTILS_H

