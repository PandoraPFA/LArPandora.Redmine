/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.h
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */
#ifndef LAR_PANDORA_HELPER_H
#define LAR_PANDORA_HELPER_H

#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"

namespace lar_pandora 
{

class LArPandoraHelper 
{
public:
    /**
     *  @brief Build a recob::Cluster object from an input vector of recob::Hit objects
     *
     *  @param id  The id code for the cluster
     *  @param hitVector  The input vector of hits
     */
    static recob::Cluster BuildCluster(const int id, const std::vector<art::Ptr<recob::Hit>> &hitVector);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_HELPER_H
