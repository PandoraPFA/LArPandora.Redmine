/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.h
 *
 *  @brief helper function for LArPandora producer modules
 *
 */
#ifndef LAR_PANDORA_HELPER_H
#define LAR_PANDORA_HELPER_H

#include "art/Persistency/Common/Ptr.h" 

#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Hit.h"

#include "Objects/ParticleFlowObject.h"

namespace lar_pandora
{

class LArPandoraHelper
{
public:
    /**
     *  @brief Build a recob::Cluster object from an input vector of recob::Hit objects
     *
     *  @param id  The id code for the cluster
     *  @param hitVector  The input vector of all hits
     *  @param hitList    The input list of isolated hits
     */
    static recob::Cluster BuildCluster(const int id, const std::vector<art::Ptr<recob::Hit>> &hitVector,
        const std::set<art::Ptr<recob::Hit>> &hitList);

    /**
     *  @brief Build a recob::Track object from an input Pandora particle flow object
     *
     *  @param id  The id code for the track
     *  @param pPfo  The particle flow object
     */
    static recob::Track BuildTrack(const int id, const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief Build a recob::Track object from an input vector of recob::SpacePoint objects
     *
     *  @param id  The id code for the track
     *  @param spacePointVector  The input vector of space points
     */
    static recob::Track BuildTrack(const int id, const std::vector<art::Ptr<recob::SpacePoint>> &spacePointVector);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_HELPER_H
