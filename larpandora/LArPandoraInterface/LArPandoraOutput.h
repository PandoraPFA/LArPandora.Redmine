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
     *  @brief  Find the index of an input object in an input list. Throw an exception if it doesn't exist
     *
     *  @param  pT the input object for which the ID should be found
     *  @param  tList a list of objects of type pT to query
     *  
     *  @return the ID of the input object
     */
    template <typename T>
    static size_t GetId(const T *const pT, const std::list<const T*> &tList);

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
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline size_t LArPandoraOutput::GetId(const T *const pT, const std::list<const T*> &tList)
{
    typename std::list<const T*>::const_iterator it(std::find(tList.begin(), tList.end(), pT));

    if (it == tList.end())
        throw cet::exception("LArPandora") << " LArPandoraOutput::GetId --- can't find the id of supplied object";

    return static_cast<size_t>(std::distance(tList.begin(), it));
}

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
