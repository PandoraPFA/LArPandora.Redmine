/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.h
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */
#ifndef LAR_PANDORA_OUTPUT_H
#define LAR_PANDORA_OUTPUT_H

#include "art/Persistency/Common/PtrMaker.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/Cluster.h"

#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlgBase.h"

#include "larpandora/LArPandoraInterface/ILArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "Pandora/PandoraInternal.h"

namespace art {class EDProducer;}
namespace pandora {class Pandora;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

class LArPandoraOutput
{
public:
    typedef std::map<size_t, std::vector<size_t> > IdToIdVectorMap;
    typedef std::map<const pandora::CaloHit *, art::Ptr<recob::Hit> > CaloHitToArtHitMap;

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

    static pandora::PfoList CollectPfos(const pandora::Pandora *const pMasterPandora);

    static void CollectPfos(const pandora::PfoList &parentPfoList, pandora::PfoList &pfoList);

    static pandora::VertexList CollectVertices(const pandora::PfoList &pfoList, IdToIdVectorMap &pfoToVerticesMap);

    static pandora::ClusterList CollectClusters(const pandora::PfoList &pfoList, IdToIdVectorMap &pfoToClustersMap);
    
    static void Collect3DHits(const pandora::ParticleFlowObject *const pPfo, pandora::CaloHitVector &caloHits);

    static pandora::CaloHitList Collect3DHits(const pandora::PfoList &pfoList, IdToIdVectorMap &pfoToThreeDHitsMap);

    static void BuildVertices(const pandora::VertexList &vertexList, std::unique_ptr< std::vector<recob::Vertex> > &outputVertices);

    static void BuildSpacePoints(const art::Event &event, const art::EDProducer *const pProducer, const pandora::CaloHitList &threeDHitList, const CaloHitToArtHitMap &caloHitToArtHitMap, std::unique_ptr< std::vector<recob::SpacePoint> > &outputSpacePoints, std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> > &outputSpacePointsToHits);

    static void BuildClusters(const art::Event &event, const art::EDProducer *const pProducer, const pandora::ClusterList &clusterList, const CaloHitToArtHitMap &pandoraHitToArtHitMap, const IdToIdVectorMap &pfoToClustersMap, std::unique_ptr< std::vector<recob::Cluster> > &outputClusters, std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > &outputClustersToHits, IdToIdVectorMap &pfoToArtClustersMap);

    static void BuildPFParticles(const art::Event &event, const art::EDProducer *const pProducer, const pandora::PfoList &pfoList, const IdToIdVectorMap &pfoToVerticesMap, const IdToIdVectorMap &pfoToThreeDHitsMap, const IdToIdVectorMap &pfoToArtClustersMap, std::unique_ptr< std::vector<recob::PFParticle> > &outputParticles, std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > &outputParticlesToVertices, std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > &outputParticlesToSpacePoints, std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> > &outputParticlesToClusters);

    static void BuildT0s(const art::Event &event, const art::EDProducer *const pProducer, const pandora::PfoList &pfoList, std::unique_ptr< std::vector<anab::T0> > &outputT0s, const CaloHitToArtHitMap &pandoraHitToArtHitMap, std::unique_ptr< art::Assns<recob::PFParticle, anab::T0> > &outputParticlesToT0s);

    static recob::PFParticle BuildPFParticle(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList);

    static recob::Vertex BuildVertex(const pandora::Vertex *const pVertex, const pandora::VertexList &vertexList);

    static void GetHitsInCluster(const pandora::Cluster *const pCluster, pandora::CaloHitVector &sortedHits);

    static std::vector<recob::Cluster> BuildClusters(const pandora::Cluster *const pCluster, const pandora::ClusterList &clusterList, const CaloHitToArtHitMap &pandoraHitToArtHitMap, IdToIdVectorMap &pandoraClusterToArtClustersMap, std::vector<HitVector> &hitVectors, size_t &nextId, cluster::ClusterParamsAlgBase &algo);

    static recob::Cluster BuildCluster(const size_t id, const HitVector &hitVector, const HitList &isolatedHits, cluster::ClusterParamsAlgBase &algo);

    static recob::SpacePoint BuildSpacePoint(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &threeDHitList);

    static bool BuildT0(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoList &pfoList, size_t &nextId, const CaloHitToArtHitMap &pandoraHitToArtHitMap, anab::T0 &t0);

    static void GetPandoraToArtHitMap(const pandora::ClusterList &clusterList, const pandora::CaloHitList &threeDHitList, const IdToHitMap &idToHitMap, CaloHitToArtHitMap &pandoraHitToArtHitMap);
    
    template <typename A, typename B>
    static void AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const size_t idA, const size_t idB, std::unique_ptr< art::Assns<A, B> > &association);

    template <typename A, typename B>
    static void AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const size_t idA, const IdToIdVectorMap &aToBMap, std::unique_ptr< art::Assns<A, B> > &association);

    template <typename A, typename B>
    static void AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const size_t idA, const std::vector< art::Ptr<B> > &bVector, std::unique_ptr< art::Assns<A, B> > &association);
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

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <typename A, typename B>
inline void LArPandoraOutput::AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const size_t idA,
    const size_t idB, std::unique_ptr< art::Assns<A, B> > &association)
{
    const art::PtrMaker<A> makePtrA(event, *pProducer);
    art::Ptr<A> pA(makePtrA(idA));

    const art::PtrMaker<B> makePtrB(event, *pProducer);
    art::Ptr<B> pB(makePtrB(idB));
    
    association->addSingle(pA, pB);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <typename A, typename B>
inline void LArPandoraOutput::AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const size_t idA,
    const IdToIdVectorMap &aToBMap, std::unique_ptr< art::Assns<A, B> > &association)
{
    IdToIdVectorMap::const_iterator it(aToBMap.find(idA));
    if (it == aToBMap.end())
        throw cet::exception("LArPandora") << " LArPandoraOutput::AddAssociation --- id doesn't exists in the assocaition map";

    const art::PtrMaker<A> makePtrA(event, *pProducer);
    art::Ptr<A> pA(makePtrA(idA));

    const art::PtrMaker<B> makePtrB(event, *pProducer);
    for (const size_t idB : it->second)
    {
        art::Ptr<B> pB(makePtrB(idB));
        association->addSingle(pA, pB);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <typename A, typename B>
inline void LArPandoraOutput::AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const size_t idA, const std::vector< art::Ptr<B> > &bVector, std::unique_ptr< art::Assns<A, B> > &association)
{
    const art::PtrMaker<A> makePtrA(event, *pProducer);
    art::Ptr<A> pA(makePtrA(idA));
    
    for (const art::Ptr<B> &pB : bVector)
        association->addSingle(pA, pB);
}

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
