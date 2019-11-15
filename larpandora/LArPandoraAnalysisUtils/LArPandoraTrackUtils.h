/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Tracks
 *
 * @author leigh.howard.whitehead@cern.ch
*/

#ifndef LAR_PANDORA_TRACK_UTILS_H
#define LAR_PANDORA_TRACK_UTILS_H

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
    static const std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Track> &track, art::Event const &evt, const std::string &label);

    /**
    * @brief Get the spacepoints associated with the track.
    *
    * @param track is the track for which we want the spacepoints
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    * 
    * @return vector of art::Ptrs to the spacepoints 
    */
    static const std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Ptr<recob::Track> &track, art::Event const &evt, const std::string &label);

    /**
    * @brief Get the particle associated with the track.
    *
    * @param track is the track for which we want the particle
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    * 
    * @return art::Ptr to the particle
    */
    static const art::Ptr<recob::PFParticle> GetParticle(const art::Ptr<recob::Track> &track, art::Event const &evt, const std::string &label);

private:

};

} // namespace lar_pandora


#endif // LAR_PANDORA_TRACK_UTILS_H

