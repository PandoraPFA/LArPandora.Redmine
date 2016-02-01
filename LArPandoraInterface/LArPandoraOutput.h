/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.h
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */
#ifndef LAR_PANDORA_OUTPUT_H
#define LAR_PANDORA_OUTPUT_H

#include "art/Persistency/Common/Ptr.h" 

#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"

#include "LArObjects/LArTrackPfo.h"

#include "LArPandoraInterface/ILArPandora.h"
#include "LArPandoraInterface/LArPandoraHelper.h"

#include <vector>
#include <set>

namespace recob {class Hit;}
namespace pandora {class Pandora; class ParticleFlowObject;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

typedef std::set< art::Ptr<recob::Hit> > HitList;

class LArPandoraOutput
{
public:
    /**
     *  @brief  Convert the Pandora PFOs into ART clusters and write into ART event
     *
     *  @param  producer the producer
     *  @param  evt the ART event
     *  @param  pPrimaryPandora the address of the primary pandora instance
     *  @param  idToHitMap the mapping from Pandora hit ID to ART hit
     *  @param  buildTracks whether to build tracks
     *  @param  buildShowers whether to build showers
     */
    static void ProduceArtOutput(art::EDProducer &producer, art::Event &evt, const pandora::Pandora *const pPrimaryPandora,
        const IdToHitMap &idToHitMap, const bool buildTracks, const bool buildShowers);

    /**
     *  @brief Build a recob::Cluster object from an input vector of recob::Hit objects
     *
     *  @param id the id code for the cluster
     *  @param hitVector the input vector of hits
     *  @param hitList the input list of isolated hits
     *  @param algo Algorithm set to fill cluster members
     *  
     *  If you don't know which algorithm to pick, StandardClusterParamsAlg is a good default.
     *  The hits that are isolated (that is, present in isolatedHits) are not fed to the cluster parameter algorithms.
     */
    static recob::Cluster BuildCluster(const int id, const HitVector &hitVector, const HitList &hitList, cluster::ClusterParamsAlgBase &algo);

    /**
     *  @brief Build a recob::Track object from an input Pandora particle flow object
     *
     *  @param id the id code for the track
     *  @param pPfo the particle flow object
     */
    static recob::Track BuildTrack(const int id, const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief Merge a set of recob::Track objects
     *
     *  @param id the id code for the tracks
     *  @param trackStateVector the vector of trajectory points for this track
     */
    static recob::Track BuildTrack(const int id, const lar_content::LArTrackStateVector &trackStateVector);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
