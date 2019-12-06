/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.h
 *
 * @brief Utility containing helpful functions for end users to access products from events
*/

#ifndef LAR_PANDORA_EVENT_UTILS_H
#define LAR_PANDORA_EVENT_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h"

#include <string>
#include <vector>

namespace lar_pandora
{
/**
 *
 * @brief LArPandoraEventUtils class
 *
*/
class LArPandoraEventUtils:LArPandoraUtilsBase
{
public:
    /**
    * @brief Get the particles from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return vector of art::Ptrs to particles
    */
    static const std::vector<art::Ptr<recob::PFParticle>> GetPFParticles(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the tracks from the event. This function shouldn't be used as the basis of an analysis.
    *
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    *
    * @return vector of art::Ptrs to tracks
    */
    static const std::vector<art::Ptr<recob::Track>> GetTracks(const art::Event &evt, const std::string &label);
    
    /**
    * @brief Get the showers from the event. This function shouldn't be used as the basis of an analysis.
    *
    * @param evt is the underlying art event
    * @param label is the label for the shower producer
    *
    * @return vector of art::Ptrs to showers
    */
    static const std::vector<art::Ptr<recob::Shower>> GetShowers(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the vertices from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the vertex producer
    *
    * @return vector of art::Ptrs to vertices
    */
    static const std::vector<art::Ptr<recob::Vertex>> GetVertices(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the spacepoints from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the spacepoint producer
    *
    * @return vector of art::Ptrs to spacepoints
    */
    static const std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the slices from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the slice producer
    *
    * @return vector of art::Ptrs to slices
    */
    static const std::vector<art::Ptr<recob::Slice>> GetSlices(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the clear cosmic ray primaries from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return vector of art::Ptrs to cosmic particles
    */
    static const std::vector<art::Ptr<recob::PFParticle>> GetClearCosmics(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the neutrino from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return atr::Ptr to the neutrino
    */
    static const art::Ptr<recob::PFParticle> GetNeutrino(const art::Event &evt, const std::string &label);

    /**
    * @brief Check to see if the event has a reconstructed neutrino
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return true if the event has a reconstructed neutrino
    */
    static const bool HasNeutrino(const art::Event &evt, const std::string &label);    
};

} // namespace lar_pandora

#endif // LAR_PANDORA_EVENT_UTILS_H

