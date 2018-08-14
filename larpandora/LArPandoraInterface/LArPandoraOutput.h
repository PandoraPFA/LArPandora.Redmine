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
    typedef std::vector<size_t> IdVector;
    typedef std::map<size_t, IdVector> IdToIdVectorMap;
    typedef std::map<const pandora::CaloHit *, art::Ptr<recob::Hit> > CaloHitToArtHitMap;
  
    typedef std::unique_ptr< std::vector<recob::PFParticle> > PFParticleCollection;
    typedef std::unique_ptr< std::vector<recob::Vertex> > VertexCollection;
    typedef std::unique_ptr< std::vector<recob::Cluster> > ClusterCollection;
    typedef std::unique_ptr< std::vector<recob::SpacePoint> > SpacePointCollection;
    typedef std::unique_ptr< std::vector<anab::T0> > T0Collection;
    typedef std::unique_ptr< std::vector<larpandoraobj::PFParticleMetadata> > PFParticleMetadataCollection;

    typedef std::unique_ptr< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> > PFParticleToMetadataCollection;
    typedef std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > PFParticleToSpacePointCollection;
    typedef std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> > PFParticleToClusterCollection;
    typedef std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > PFParticleToVertexCollection;
    typedef std::unique_ptr< art::Assns<recob::PFParticle, anab::T0> > PFParticleToT0Collection;

    typedef std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > ClusterToHitCollection;
    typedef std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> > SpacePointToHitCollection;

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

        /**
         *  @brief  Check the parameters and throw an exception if they are not valid
         */
        void Validate() const;

        const pandora::Pandora *m_pPrimaryPandora;              ///<
        art::EDProducer        *m_pProducer;                    ///<
        bool                    m_shouldRunStitching;           ///<
        bool                    m_shouldProduceAllOutcomes;     ///< If all outcomes should be produced in separate collections (choose false if you only require the consolidated output)
        std::string             m_allOutcomesInstanceLabel;     ///< The label for the instance producing all outcomes
    };

    /**
     *  @brief  Convert the Pandora PFOs into ART clusters and write into ART event
     *
     *  @param  settings the settings
     *  @param  idToHitMap the mapping from Pandora hit ID to ART hit
     *  @param  evt the ART event
     */
    static void ProduceArtOutput(const Settings &settings, const IdToHitMap &idToHitMap, art::Event &evt);

private: 
    /**
     *  @brief  Collect the current pfos (including all downstream pfos) from the master pandora instance
     *
     *  @param  pPrimaryPandora address of master pandora instance
     *
     *  @return a sorted list of all pfos to convert to ART PFParticles
     */
    static pandora::PfoVector CollectPfos(const pandora::Pandora *const pPrimaryPandora);
    
    /**
     *  @brief  Collect the pfos (including all downstream pfos) from the master and daughter pandora instances
     *
     *  @param  pPrimaryPandora address of master pandora instance
     *
     *  @return a sorted list of all pfos to convert to ART PFParticles
     */
    static pandora::PfoVector CollectAllPfoOutcomes(const pandora::Pandora *const pPrimaryPandora);

    /**
     *  @brief  Collect a sorted list of all downstream pfos of an input list of parent
     *
     *  @param  parentPfoList the input list of parent pfos
     *  @param  pfoList the sorted output list of all downstream pfos
     */
    static void CollectPfos(const pandora::PfoList &parentPfoList, pandora::PfoVector &pfoVector);

    /**
     *  @brief  Collect all vertices contained in the input pfo list
     *          Order is guaranteed provided pfoList is ordered
     *
     *  @param  pfoList the input list of pfos
     *  @param  pfoToVerticesMap the output mapping from pfo ID to vertex IDs (zero or one)
     *
     *  @return the list of vertices collected
     */
    static pandora::VertexVector CollectVertices(const pandora::PfoVector &pfoVector, IdToIdVectorMap &pfoToVerticesMap);

    /**
     *  @brief  Collect a sorted list of all 2D clusters contained in the input pfo list
     *          Order is guaranteed provided pfoList is ordered
     *
     *  @param  pfoList the input list of pfos
     *  @param  pfoToClustersMap the output mapping from pfo ID to cluster IDs
     *
     *  @return the list of clusters collected
     */
    static pandora::ClusterList CollectClusters(const pandora::PfoVector &pfoVector, IdToIdVectorMap &pfoToClustersMap);
    
    /**
     *  @brief  Collect a sorted vector of all 3D hits in the input pfo
     *
     *  @param  pPfo the input pfo
     *  @param  caloHits the sorted output vector of 3D hits
     */
    static void Collect3DHits(const pandora::ParticleFlowObject *const pPfo, pandora::CaloHitVector &caloHits);

    /**
     *  @brief  Collect a sorted list of all 3D hits contained in the input pfo list
     *          Order is guaranteed provided pfoList is ordered
     *
     *  @param  pfoList the input list of pfos
     *  @param  pfoToThreeDHitsMap the output mapping from pfo ID to 3D hit IDs
     *
     *  @return the list of 3D hits collected
     */
    static pandora::CaloHitList Collect3DHits(const pandora::PfoVector &pfoVector, IdToIdVectorMap &pfoToThreeDHitsMap);

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
    
    template <typename T>
    static size_t GetId(const T *const pT, const std::vector<const T*> &tVector);
    
    /**
     *  @brief  Collect all 2D and 3D hits that were used / produced in the reconstruction and map them to their corresponding ART hit
     *
     *  @param  clusterList input list of all 2D clusters to be output
     *  @param  threeDHitList input list of all 3D hits to be output (as spacepoints)
     *  @param  idToHitMap input mapping from pandora hit ID to ART hit
     *  @param  pandoraHitToArtHitMap output mapping from pandora hit to ART hit
     */
    static void GetPandoraToArtHitMap(const pandora::ClusterList &clusterList, const pandora::CaloHitList &threeDHitList,
        const IdToHitMap &idToHitMap, CaloHitToArtHitMap &pandoraHitToArtHitMap);

    /**
     *  @brief  Look up ART hit from an input Pandora hit
     *
     *  @param  idToHitMap the mapping between Pandora and ART hits
     *  @param  pCaloHit the input Pandora hit (2D)
     */
    static art::Ptr<recob::Hit> GetHit(const IdToHitMap &idToHitMap, const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief  Convert pandora vertices to ART vertices and add them to the output vector
     *
     *  @param  vertexList the input list of pandora vertices
     *  @param  outputVertices the output vector of ART vertices
     */
    static void BuildVertices(const pandora::VertexVector &vertexVector, VertexCollection &outputVertices);

    /**
     *  @brief  Convert pandora 3D hits to ART spacepoints and add them to the output vector
     *          Create the associations between spacepoints and hits
     *
     *  @param  event the art event
     *  @param  pProducer the address of the producer module
     *  @param  threeDHitList the input list of 3D hits to convert
     *  @param  pandoraHitToArtHitMap the input mapping from pandora hits to ART hits
     *  @param  outputSpacePoints the output vector of spacepoints
     *  @param  outputSpacePointsToHits the output associations between spacepoints and hits
     */
    static void BuildSpacePoints(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const pandora::CaloHitList &threeDHitList, const CaloHitToArtHitMap &pandoraHitToArtHitMap, SpacePointCollection &outputSpacePoints,
        SpacePointToHitCollection &outputSpacePointsToHits);

    /**
     *  @brief  Convert pandora 2D clusters to ART clusters and add them to the output vector
     *          Create the associations between clusters and hits.
     *          For multiple drift volumes, each pandora cluster can correspond to multiple ART clusters.
     *
     *  @param  event the art event
     *  @param  pProducer the address of the producer module
     *  @param  clusterList the input list of 2D pandora clusters to convert
     *  @param  pandoraHitToArtHitMap the input mapping from pandora hits to ART hits
     *  @param  pfoToClustersMap the input mapping from pfo ID to cluster IDs
     *  @param  outputClusters the output vector of clusters
     *  @param  outputClustersToHits the output associations between clusters and hits
     *  @param  pfoToArtClustersMap the output mapping from pfo ID to art cluster ID
     */
    static void BuildClusters(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const pandora::ClusterList &clusterList, const CaloHitToArtHitMap &pandoraHitToArtHitMap, const IdToIdVectorMap &pfoToClustersMap,
        ClusterCollection &outputClusters, ClusterToHitCollection &outputClustersToHits, IdToIdVectorMap &pfoToArtClustersMap);

    /**
     *  @brief  Convert between pfos and PFParticles and add them to the output vector
     *          Create the associations between PFParticle and vertices, spacepoints and clusters
     *
     *  @param  event the art event
     *  @param  pProducer the address of the producer module
     *  @param  pfoList the input list of pfos to convert
     *  @param  pfoToVerticesMap the input mapping from pfo ID to vertex IDs
     *  @param  pfoToThreeDHitsMap the input mapping from pfo ID to 3D hit IDs
     *  @param  pfoToArtClustersMap the input mapping from pfo ID to ART cluster IDs
     *  @param  outputParticle the output vector of PFParticles
     *  @param  outputParticlesToVertices the output associations between PFParticles and vertices
     *  @param  outputParticlesToSpacePoints the output associations between PFParticles and spacepoints
     *  @param  outputParticlesToClusters the output associations between PFParticles and clusters
     */
    static void BuildPFParticles(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const pandora::PfoVector &pfoVector, const IdToIdVectorMap &pfoToVerticesMap, const IdToIdVectorMap &pfoToThreeDHitsMap, 
        const IdToIdVectorMap &pfoToArtClustersMap, PFParticleCollection &outputParticles, 
        PFParticleToVertexCollection &outputParticlesToVertices, PFParticleToSpacePointCollection &outputParticlesToSpacePoints,
        PFParticleToClusterCollection &outputParticlesToClusters);

    /**
     *  @brief  Build metadata objects from a list of input pfos
     *
     *  @param  event the art event
     *  @param  pProducer the address of the producer module
     *  @param  pfoList the input list of pfos
     *  @param  outputParticleMetadata the output vector of PFParticleMetadata
     *  @param  outputParticlesToMetadata the output associations between PFParticles and metadata
     */
    static void BuildParticleMetadata(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel, 
        const pandora::PfoVector &pfoVector, PFParticleMetadataCollection &outputParticleMetadata,
        PFParticleToMetadataCollection &outputParticlesToMetadata);

    /**
     *  @brief  Calculate the T0 of each pfos and add them to the output vector
     *          Create the associations between PFParticle and T0s
     *
     *  @param  event the art event
     *  @param  pProducer the address of the producer module
     *  @param  pfoList the input list of pfos
     *  @param  outputT0s the output vector of T0s
     *  @param  pandoraHitToArtHitMap the input mapping from pandora hits to ART hits
     *  @param  outputParticlesToT0s the output associations between PFParticles and T0s
     */
    static void BuildT0s(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const pandora::PfoVector &pfoVector, T0Collection &outputT0s, const CaloHitToArtHitMap &pandoraHitToArtHitMap,
        PFParticleToT0Collection &outputParticlesToT0s);

    /**
     *  @brief  Convert from a pandora vertex to an ART vertex
     *
     *  @param  pVertex the input vertex
     *  @param  vertexId the id of the vertex to produce
     *
     *  @param  the ART vertex
     */
    static recob::Vertex BuildVertex(const pandora::Vertex *const pVertex, const size_t vertexId);
    
    /**
     *  @brief  Convert from a pandora 3D hit to an ART spacepoint
     *
     *  @param  pCaloHit the input hit
     *  @param  spacePointId the id of the space-point to produce
     *
     *  @param  the ART spacepoint
     */
    static recob::SpacePoint BuildSpacePoint(const pandora::CaloHit *const pCaloHit, const size_t spacePointId);

    /**
     *  @brief  Collect a sorted list of all 2D hits in a cluster
     *
     *  @param  pCluster the input cluster
     *  @param  sortedHits the output vector of sorted 2D hits
     */
    static void GetHitsInCluster(const pandora::Cluster *const pCluster, pandora::CaloHitVector &sortedHits);

    /**
     *  @brief  Convert from a pandora 2D cluster to a vector of ART clusters (produce multiple if the cluster is split over drift volumes)
     *
     *  @param  pCluster the input cluster
     *  @param  clusterList the input list of clusters
     *  @param  pandoraHitToArtHitMap the input mapping from pandora hits to ART hits
     *  @param  pandoraClusterToArtClustersMap output mapping from pandora cluster ID to art cluster IDs
     *  @param  hitVectors the output vectors of hits for each cluster produced used to produce associations
     *  @param  algo algorithm set to fill cluster members
     *
     *  @param  the vector of ART clusters
     */
    static std::vector<recob::Cluster> BuildClusters(const pandora::Cluster *const pCluster, const pandora::ClusterList &clusterList,
        const CaloHitToArtHitMap &pandoraHitToArtHitMap, IdToIdVectorMap &pandoraClusterToArtClustersMap, std::vector<HitVector> &hitVectors,
        size_t &nextId, cluster::ClusterParamsAlgBase &algo);

    /**
     *  @brief  Build an ART cluster from an input vector of ART hits
     *
     *  @param  id the id code for the cluster
     *  @param  hitVector the input vector of hits
     *  @param  isolatedHits the input list of isolated hits
     *  @param  algo algorithm set to fill cluster members
     * 
     *  @return the ART cluster
     *
     *  If you don't know which algorithm to pick, StandardClusterParamsAlg is a good default.
     *  The hits that are isolated (that is, present in isolatedHits) are not fed to the cluster parameter algorithms.
     */
    static recob::Cluster BuildCluster(const size_t id, const HitVector &hitVector, const HitList &isolatedHits,
        cluster::ClusterParamsAlgBase &algo);
    
    /**
     *  @brief  Convert from a pfo to and ART PFParticle
     *
     *  @param  pPfo the input pfo to convert
     *  @param  pfoId the id of the pfo to produce
     *  @param  pfoList the input list of pfos
     *
     *  @param  the ART PFParticle
     */
    static recob::PFParticle BuildPFParticle(const pandora::ParticleFlowObject *const pPfo, const size_t pfoId, const pandora::PfoVector &pfoVector);

    /**
     *  @brief  If required, build a T0 for the input pfo
     *
     *  @param  pPfo the input pfo
     *  @param  pfoList the input list of pfos
     *  @param  nextId the ID of the T0 - will be incremented if the t0 was produced
     *  @param  pandoraHitToArtHitMap the input mapping from pandora hits to ART hits
     *  @param  t0 the output T0
     *
     *  @return if a T0 was produced (calculated from the stitching hit shift distance)
     */
    static bool BuildT0(const pandora::ParticleFlowObject *const pPfo, const pandora::PfoVector &pfoVector, size_t &nextId,
        const CaloHitToArtHitMap &pandoraHitToArtHitMap, anab::T0 &t0);
    
    /**
     *  @brief  Convert X0 correction into T0 correction
     *
     *  @param  hit the input ART hit
     *  @param  pCaloHit the output Pandora hit
     *
     *  @return T0 relative to input hit in nanoseconds
     */
    static double CalculateT0(const art::Ptr<recob::Hit> hit, const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief  Add an association between objects with two given ids
     *
     *  @param  event the ART event
     *  @param  pProducer the address of the producer module
     *  @param  idA the id of an object of type A
     *  @param  idB the id of an object of type B to associate to the first object
     *  @param  association the output association to update
     */
    template <typename A, typename B>
    static void AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const size_t idA, const size_t idB, std::unique_ptr< art::Assns<A, B> > &association);

    /**
     *  @brief  Add associations between input objects
     *
     *  @param  event the ART event
     *  @param  pProducer the address of the producer module
     *  @param  idA the id of an object of type A
     *  @param  aToBMap the input mapping from IDs of objects of type A to IDs of objects of type B to associate
     *  @param  association the output association to update
     */
    template <typename A, typename B>
    static void AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const size_t idA, const IdToIdVectorMap &aToBMap, std::unique_ptr< art::Assns<A, B> > &association);

    /**
     *  @brief  Add associations between input objects
     *
     *  @param  event the ART event
     *  @param  pProducer the address of the producer module
     *  @param  idA the id of an object of type A
     *  @param  bVector the input vector of IDs of objects of type B to associate
     *  @param  association the output association to update
     */
    template <typename A, typename B>
    static void AddAssociation(const art::Event &event, const art::EDProducer *const pProducer, const std::string &instanceLabel,
        const size_t idA, const std::vector< art::Ptr<B> > &bVector, std::unique_ptr< art::Assns<A, B> > &association);
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

template <typename T>
inline size_t LArPandoraOutput::GetId(const T *const pT, const std::vector<const T*> &tVector)
{
    typename std::vector<const T*>::const_iterator it(std::find(tVector.begin(), tVector.end(), pT));

    if (it == tVector.end())
        throw cet::exception("LArPandora") << " LArPandoraOutput::GetId --- can't find the id of supplied object";

    return static_cast<size_t>(std::distance(tVector.begin(), it));
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <typename A, typename B>
inline void LArPandoraOutput::AddAssociation(const art::Event &event, const art::EDProducer *const pProducer,
    const std::string &instanceLabel, const size_t idA, const size_t idB, std::unique_ptr< art::Assns<A, B> > &association)
{
    const art::PtrMaker<A> makePtrA(event, *pProducer, instanceLabel);
    art::Ptr<A> pA(makePtrA(idA));

    const art::PtrMaker<B> makePtrB(event, *pProducer, instanceLabel);
    art::Ptr<B> pB(makePtrB(idB));
    
    association->addSingle(pA, pB);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <typename A, typename B>
inline void LArPandoraOutput::AddAssociation(const art::Event &event, const art::EDProducer *const pProducer,
    const std::string &instanceLabel, const size_t idA, const IdToIdVectorMap &aToBMap, std::unique_ptr< art::Assns<A, B> > &association)
{
    IdToIdVectorMap::const_iterator it(aToBMap.find(idA));
    if (it == aToBMap.end())
        throw cet::exception("LArPandora") << " LArPandoraOutput::AddAssociation --- id doesn't exists in the assocaition map";

    const art::PtrMaker<A> makePtrA(event, *pProducer, instanceLabel);
    art::Ptr<A> pA(makePtrA(idA));

    const art::PtrMaker<B> makePtrB(event, *pProducer, instanceLabel);
    for (const size_t idB : it->second)
    {
        art::Ptr<B> pB(makePtrB(idB));
        association->addSingle(pA, pB);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <typename A, typename B>
inline void LArPandoraOutput::AddAssociation(const art::Event &event, const art::EDProducer *const pProducer,
    const std::string &instanceLabel, const size_t idA, const std::vector< art::Ptr<B> > &bVector,
    std::unique_ptr< art::Assns<A, B> > &association)
{
    const art::PtrMaker<A> makePtrA(event, *pProducer, instanceLabel);
    art::Ptr<A> pA(makePtrA(idA));
    
    for (const art::Ptr<B> &pB : bVector)
        association->addSingle(pA, pB);
}

} // namespace lar_pandora

#endif //  LAR_PANDORA_OUTPUT_H
