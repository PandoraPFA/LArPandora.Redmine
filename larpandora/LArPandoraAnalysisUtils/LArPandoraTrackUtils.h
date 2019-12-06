/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Tracks
*/

#ifndef LAR_PANDORA_TRACK_UTILS_H
#define LAR_PANDORA_TRACK_UTILS_H

#include "art/Framework/Principal/Event.h"

// Access the type defs defined in the helper
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h"

#include <string>
#include <vector>

namespace anab
{
class Calorimetry;
}

namespace lar_pandora
{
/**
 *
 * @brief LArPandoraTrackUtils class
 *
*/
class LArPandoraTrackUtils:LArPandoraUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the track.
    *
    * @param track is the track for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    * 
    * @return vector of art::Ptrs to the hits 
    */
    static const std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the spacepoints associated with the track.
    *
    * @param track is the track for which we want the spacepoints
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    * 
    * @return vector of art::Ptrs to the spacepoints 
    */
    static const std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the particle associated with the track.
    *
    * @param track is the track for which we want the particle
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    * 
    * @return art::Ptr to the particle
    */
    static const art::Ptr<recob::PFParticle> GetPFParticle(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the particle associated with the track.
    *
    * @param track is the track for which we want the particle
    * @param evt is the underlying art event
    * @param trackLabel is the label for the track producer
    * @param caloLabel is the label for the calorimetry producer
    *
    * @return art::Ptr to the particle
    */
    static const art::Ptr<anab::Calorimetry> GetCalorimetry(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &trackLabel, const std::string &caloLabel);
};

} // namespace lar_pandora

#endif // LAR_PANDORA_TRACK_UTILS_H

