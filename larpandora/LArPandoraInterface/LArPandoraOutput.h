/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.h
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */
#ifndef LAR_PANDORA_OUTPUT_H
#define LAR_PANDORA_OUTPUT_H

#include "lardataobj/RecoBase/Cluster.h"

#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"

#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

#include "larpandora/LArPandoraInterface/ILArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace art {class EDProducer;}
namespace pandora {class Pandora; class CaloHit;}

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
        bool                    m_shouldRunStitching;           ///<
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

    /**
     *  @brief Convert X0 correction into T0 correction
     *
     *  @param hit the input ART hit
     *  @param pCaloHit the output Pandora hit
     *
     *  @return T0 relative to input hit in nanoseconds
     */
    static double CalculateT0(const art::Ptr<recob::Hit> hit, const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief Get metadata about a given pfo
     *
     *  @param pPfo the Pandora pfo
     *
     *  @return metadata collected about the input pfo
     */
    static larpandoraobj::PFParticleMetadata GetPFParticleMetadata(const pandora::ParticleFlowObject *const pPfo);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
