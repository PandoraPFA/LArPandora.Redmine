/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.h
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */
#ifndef LAR_PANDORA_OUTPUT_H
#define LAR_PANDORA_OUTPUT_H

#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"

#include "LArObjects/LArTrackPfo.h"

#include "larpandora/LArPandoraInterface/ILArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace art {class EDProducer;}
namespace pandora {class Pandora; class ParticleFlowObject;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

class LArPandoraOutput
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        const pandora::Pandora *m_pPrimaryPandora;              ///< 
        art::EDProducer        *m_pProducer;                    ///< 
        bool                    m_buildTracks;                  ///<
        bool                    m_buildShowers;                 ///<
        bool                    m_buildStitchedParticles;       ///<
        bool                    m_buildSingleVolumeParticles;   ///<
    };

    /**
     *  @brief  Convert the Pandora PFOs into ART clusters and write into ART event
     *
     *  @param  settings the settings
     *  @param  idToHitMap the mapping from Pandora hit ID to ART hit
     *  @param  evt the ART event
     */
    static void ProduceArtOutput(const Settings &settings, const IdToHitMap &idToHitMap, art::Event &evt);

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
     *  @param pTrackStateVector the vector of trajectory points for this track
     */
    static recob::Track BuildTrack(const int id, const lar_content::LArTrackStateVector *const pTrackStateVector);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
