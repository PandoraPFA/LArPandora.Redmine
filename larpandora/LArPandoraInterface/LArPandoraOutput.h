/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.h
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */
#ifndef LAR_PANDORA_OUTPUT_H
#define LAR_PANDORA_OUTPUT_H

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"

#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"
#include "larreco/Calorimetry/LinearEnergyAlg.h"

#include "larpandoracontent/LArObjects/LArTrackPfo.h" // For track state definitions
#include "larpandoracontent/LArObjects/LArShowerPfo.h" // For shower parameters

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
        unsigned int            m_minTrajectoryPoints;          ///<
        bool                    m_buildShowers;                 ///<
        bool                    m_buildStitchedParticles;       ///<
        bool                    m_buildSingleVolumeParticles;   ///<
        calo::LinearEnergyAlg const* m_showerEnergyAlg;         ///<
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
     *  @param isolatedHits the input list of isolated hits
     *  @param algo Algorithm set to fill cluster members
     *  
     *  If you don't know which algorithm to pick, StandardClusterParamsAlg is a good default.
     *  The hits that are isolated (that is, present in isolatedHits) are not fed to the cluster parameter algorithms.
     */
    static recob::Cluster BuildCluster(const int id, const HitVector &hitVector, const HitList &isolatedHits, cluster::ClusterParamsAlgBase &algo);

    /**
     *  @brief Check (valid) trajectory points is at least the minimum
     *
     *  @param pTrackStateVector the vector of trajectory points for this track
     */
    static bool MinTrajectoryPoints(const lar_content::LArTrackStateVector *const pTrackStateVector, const unsigned int minTrajectoryPoints);

    /**
     *  @brief Build a recob::Track object
     *
     *  @param id the id code for the track
     *  @param pTrackStateVector the vector of trajectory points for this track
     */
    static recob::Track BuildTrack(const int id, const lar_content::LArTrackStateVector *const pTrackStateVector);

    /**
     *  @brief Build a recob::Shower object
     *
     *  @param pLArShowerPfo the object of the shower parameters filled in pandora
     */
    static recob::Shower BuildShower(const lar_content::LArShowerPfo *const pLArShowerPfo);

    /**
     *  @brief Build a recob::PCAxis object
     *
     *  @param pLArShowerPfo the object of the shower parameters filled in pandora
     */
    static recob::PCAxis BuildShowerPCA(const lar_content::LArShowerPfo *const pLArShowerPfo);

    /**
     *  @brief Build a reco::Seed object
     *
     *  @param trackState the trajectory point for this seed
     */
    static recob::Seed BuildSeed(const lar_content::LArTrackState &trackState);

    /**
     *  @brief Build a recob::SpacePoint object
     *
     *  @param id the id code for the spacepoint
     *  @param pCaloHit the input Pandora hit (3D)
     */
    static recob::SpacePoint BuildSpacePoint(const int id, const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief Lookup ART hit from an input Pandora hit
     *
     *  @param idToHitMap the mapping between Pandora and ART hits
     *  @param pCaloHit the input Pandora hit (2D)
     */
    static art::Ptr<recob::Hit> GetHit(const IdToHitMap &idToHitMap, const pandora::CaloHit *const pCaloHit);

};

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
