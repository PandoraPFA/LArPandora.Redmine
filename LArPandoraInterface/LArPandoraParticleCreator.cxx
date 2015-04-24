// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "Utilities/AssociationUtil.h"
#include "SimulationBase/MCTruth.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/PFParticle.h"
#include "RecoBase/Vertex.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"

// Pandora includes
#include "Objects/ParticleFlowObject.h"
#include "Objects/Vertex.h"
#include "LArContent.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

// Local LArPandora includes
#include "LArPandoraInterface/LArPandoraParticleCreator.h"
#include "LArPandoraInterface/LArPandoraHelper.h"

// System includes
#include <iostream>
#include <limits>

namespace lar_pandora {

LArPandoraParticleCreator::LArPandoraParticleCreator(fhicl::ParameterSet const &pset) : LArPandoraBase(pset)
{
    m_buildTracks = pset.get<bool>("BuildTracks", false);
    m_buildShowers = pset.get<bool>("BuildShowers", false);

    produces< std::vector<recob::PFParticle> >();
    produces< std::vector<recob::SpacePoint> >();
    produces< std::vector<recob::Cluster> >(); 
    produces< std::vector<recob::Seed> >();
    produces< std::vector<recob::Vertex> >();

    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();

    if (m_buildTracks)
    {
        produces< std::vector<recob::Track> >(); 
        produces< art::Assns<recob::PFParticle, recob::Track> >();
        produces< art::Assns<recob::Track, recob::Hit> >();
    }

    if (m_buildShowers)
    {
        produces< std::vector<recob::Shower> >(); 
        produces< art::Assns<recob::PFParticle, recob::Shower> >();
        produces< art::Assns<recob::Shower, recob::Hit> >();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraParticleCreator::~LArPandoraParticleCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraParticleCreator::beginJob()
{
    if (m_enableMonitoring)
        this->InitializeMonitoring();

    return LArPandoraBase::beginJob();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraParticleCreator::produce(art::Event &evt)
{ 
    mf::LogInfo("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] *** " << std::endl;

    cet::cpu_timer theClock;

    if (m_enableMonitoring)
    {   
        theClock.start();
    }

    // Collect ART objects
    HitVector artHits;
    SimChannelVector artSimChannels;
    HitsToTrackIDEs artHitsToTrackIDEs;
    MCParticleVector artMCParticleVector;
    MCTruthToMCParticles artMCTruthToMCParticles;
    MCParticlesToMCTruth artMCParticlesToMCTruth;

    LArPandoraCollector::CollectHits(evt, m_hitfinderModuleLabel, artHits);

    if (m_enableMCParticles && !evt.isRealData())
    {
        LArPandoraCollector::CollectMCParticles(evt, m_geantModuleLabel, artMCParticleVector);
        LArPandoraCollector::CollectMCParticles(evt, m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth);
        LArPandoraCollector::CollectSimChannels(evt, m_geantModuleLabel, artSimChannels);
        LArPandoraCollector::BuildMCParticleHitMaps(artHits, artSimChannels, artHitsToTrackIDEs);
    }

    if (m_enableMonitoring)
    { 
        theClock.stop();
        m_collectionTime = theClock.accumulated_real_time();
        theClock.reset();
        theClock.start();
    }

    // Create PANDORA objects
    HitMap pandoraHits;

    this->CreatePandoraHits2D(artHits, pandoraHits);

    if (m_enableMCParticles && !evt.isRealData())
    {
        this->CreatePandoraParticles(artMCTruthToMCParticles, artMCParticlesToMCTruth);
        this->CreatePandoraParticles2D(artMCParticleVector);
        this->CreatePandoraLinks2D(pandoraHits, artHitsToTrackIDEs);
    }

    if (m_enableMonitoring)
    { 
        theClock.stop();
        m_inputTime = theClock.accumulated_real_time();
        theClock.reset();
        theClock.start();
    }

    // Run PANDORA algorithms
    this->RunPandoraInstances();

    if (m_enableMonitoring)
    { 
        theClock.stop();
        m_processTime = theClock.accumulated_real_time(); 
        theClock.reset();
        theClock.start();
    }

    // Output ART objects and Reset PANDORA algorithms
    if (m_enableProduction)
    {
        this->ProduceArtOutput(evt, pandoraHits);
    }

    this->ResetPandoraInstances();

    if (m_enableMonitoring)
    {   theClock.stop();
        m_outputTime = theClock.accumulated_real_time(); 
        theClock.reset();
        theClock.start();
    }

    // Write out monitoring information
    if (m_enableMonitoring)
    {
        m_run   = evt.run();
        m_event = evt.id().event();
        m_hits = static_cast<int>(artHits.size());
        m_pRecoTree->Fill();
    }

    mf::LogDebug("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "]  Done! *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraParticleCreator::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("hits", &m_hits, "hits/I");
    m_pRecoTree->Branch("collectionTime", &m_collectionTime, "collectionTime/F");
    m_pRecoTree->Branch("inputTime", &m_inputTime, "inputTime/F");
    m_pRecoTree->Branch("processTime", &m_processTime, "processTime/F");
    m_pRecoTree->Branch("outputTime", &m_outputTime, "outputTime/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraParticleCreator::ProduceArtOutput(art::Event &evt, const HitMap &hitMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ProduceArtOutput() *** " << std::endl;

    // Concatenate output PANDORA particles
    pandora::PfoList concatenatedPfoList;

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); pIter != pIterEnd;
        ++pIter)
    {
        const pandora::Pandora *const pPandora = pIter->second;

        // Get list of Pandora particles (ATTN: assume that all reco particles live in curent list)
        const pandora::PfoList *pPfoList = NULL;
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pPandora, pPfoList));

        if (NULL == pPfoList)
            continue;

        concatenatedPfoList.insert(pPfoList->begin(), pPfoList->end());
    }

    if (concatenatedPfoList.empty())
        mf::LogDebug("LArPandora") << "   Warning: No reconstructed particles for this event " << std::endl;


    // Set up ART outputs
    std::unique_ptr< std::vector<recob::PFParticle> > outputParticles( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::SpacePoint> > outputSpacePoints( new std::vector<recob::SpacePoint> );
    std::unique_ptr< std::vector<recob::Cluster> >    outputClusters( new std::vector<recob::Cluster> );
    std::unique_ptr< std::vector<recob::Seed> >       outputSeeds( new std::vector<recob::Seed> );
    std::unique_ptr< std::vector<recob::Vertex> >     outputVertices( new std::vector<recob::Vertex> );
    std::unique_ptr< std::vector<recob::Track> >      outputTracks( new std::vector<recob::Track> );
    std::unique_ptr< std::vector<recob::Shower> >     outputShowers( new std::vector<recob::Shower> );

    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >    outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> >       outputParticlesToSeeds( new art::Assns<recob::PFParticle, recob::Seed> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> >     outputParticlesToVertices( new art::Assns<recob::PFParticle, recob::Vertex> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> >      outputParticlesToTracks( new art::Assns<recob::PFParticle, recob::Track> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower> >     outputParticlesToShowers( new art::Assns<recob::PFParticle, recob::Shower> );
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> >             outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit> >            outputShowersToHits( new art::Assns<recob::Shower, recob::Hit> );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> >        outputSpacePointsToHits( new art::Assns<recob::SpacePoint, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> >           outputClustersToHits( new art::Assns<recob::Cluster, recob::Hit> );

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    cluster::StandardClusterParamsAlg ClusterParamAlgo;
    

    // Obtain a sorted vector of all output Pfos and their daughters
    pandora::PfoList connectedPfoList;
    lar_content::LArPfoHelper::GetAllConnectedPfos(concatenatedPfoList, connectedPfoList);
    pandora::PfoVector pfoVector(connectedPfoList.begin(), connectedPfoList.end());
    std::sort(pfoVector.begin(), pfoVector.end(), lar_content::LArPfoHelper::SortByNHits);

    typedef std::map< const pandora::ParticleFlowObject*, size_t> ThreeDParticleMap;
    typedef std::map< const pandora::Vertex*, unsigned int> ThreeDVertexMap; 
    typedef std::set< art::Ptr<recob::Hit> > HitList;
    typedef std::map< int, HitVector > HitArray;

    int vertexCounter(0);
    int spacePointCounter(0);
    int clusterCounter(0);
    int trackCounter(0);
    size_t particleCounter(0);

    // Build maps of Pandora particles and Pandora vertices
    pandora::VertexVector vertexVector;
    ThreeDParticleMap particleMap;
    ThreeDVertexMap vertexMap;

    for (pandora::PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const pandora::ParticleFlowObject *const pPfo = *pIter;
        particleMap.insert( std::pair<const pandora::ParticleFlowObject*, size_t>(pPfo, particleCounter++) ); 

        if (!pPfo->GetVertexList().empty())
        {
            if(pPfo->GetVertexList().size() != 1)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const pandora::Vertex *const pVertex = *(pPfo->GetVertexList().begin());

            if (vertexMap.end() != vertexMap.find(pVertex))
                continue;

            vertexVector.push_back(pVertex);
            vertexMap.insert( std::pair<const pandora::Vertex*, unsigned int>(pVertex, vertexVector.size() - 1) ); 
        }
    }

    // Loop over Pandora vertices and build recob::Vertices
    for (pandora::VertexVector::iterator vIter = vertexVector.begin(), vIterEnd = vertexVector.end(); vIter != vIterEnd; ++vIter)
    {
        const pandora::Vertex *const pVertex = *vIter;

        ThreeDVertexMap::const_iterator wIter = vertexMap.find(pVertex);
        if (vertexMap.end() == wIter)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const pandora::CartesianVector vtxPos(pVertex->GetPosition());
        double pos[3] = { vtxPos.GetX(), vtxPos.GetY(), vtxPos.GetZ() };
        
        recob::Vertex newVertex(pos, vertexCounter++);
        outputVertices->push_back(newVertex);
    }

    // Loop over Pandora particles and build recob::PFParticles
    for (pandora::PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const pandora::ParticleFlowObject *const pPfo = *pIter;

        // Get Pfo ID
        ThreeDParticleMap::const_iterator qIter = particleMap.find(pPfo);
        if (particleMap.end() == qIter)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const size_t pfoIdCode(qIter->second);

        // Get Pfo Parents
        size_t parentIdCode(recob::PFParticle::kPFParticlePrimary);
        const pandora::PfoList &parentList(pPfo->GetParentPfoList());

        if (!parentList.empty())
        {
            if (parentList.size() != 1)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            for (pandora::PfoList::const_iterator parentIter = parentList.begin(), parentIterEnd = parentList.end();
                parentIter != parentIterEnd; ++parentIter)
            {
                const pandora::ParticleFlowObject *const pParentPfo = *parentIter;

                ThreeDParticleMap::const_iterator parentIdIter = particleMap.find(pParentPfo);
                if (particleMap.end() == parentIdIter)
                    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

                parentIdCode = parentIdIter->second;
            }
        }

        // Get Pfo Daughters
        std::vector<size_t> daughterIdCodes;
        const pandora::PfoList &daughterList(pPfo->GetDaughterPfoList());

        if (!daughterList.empty())
        {
            for (pandora::PfoList::const_iterator daughterIter = daughterList.begin(), daughterIterEnd = daughterList.end();
                daughterIter != daughterIterEnd; ++daughterIter)
            {
                const pandora::ParticleFlowObject *const pDaughterPfo = *daughterIter;

                ThreeDParticleMap::const_iterator daughterIdIter = particleMap.find(pDaughterPfo);
                if (particleMap.end() == daughterIdIter)
                    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

                const size_t daughterIdCode(daughterIdIter->second);
                daughterIdCodes.push_back(daughterIdCode);
            }
        }

        // Build Particle
        recob::PFParticle newParticle(pPfo->GetParticleId(), pfoIdCode, parentIdCode, daughterIdCodes);
        outputParticles->push_back(newParticle);

        // Build 3D Space Points 
        pandora::CaloHitList pandoraHitList3D;
        lar_content::LArPfoHelper::GetCaloHits(pPfo, pandora::TPC_3D, pandoraHitList3D);
        pandora::CaloHitVector pandoraHitVector3D(pandoraHitList3D.begin(), pandoraHitList3D.end());
        std::sort(pandoraHitVector3D.begin(), pandoraHitVector3D.end(), lar_content::LArClusterHelper::SortByPosition);

        for (pandora::CaloHitVector::const_iterator hIter = pandoraHitVector3D.begin(), hIterEnd = pandoraHitVector3D.end(); hIter != hIterEnd; ++hIter)
        {
            const pandora::CaloHit *const pCaloHit3D = *hIter;

            if (pandora::TPC_3D != pCaloHit3D->GetHitType())
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const pandora::CaloHit *const pCaloHit2D = static_cast<const pandora::CaloHit*>(pCaloHit3D->GetParentCaloHitAddress());
            const void *const pHitAddress(pCaloHit2D->GetParentCaloHitAddress());
            const intptr_t hitID_temp((intptr_t)(pHitAddress)); // TODO
            const int hitID((int)(hitID_temp));

            HitMap::const_iterator artIter = hitMap.find(hitID);

            if (hitMap.end() == artIter)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            HitVector hitVector;
            const art::Ptr<recob::Hit> hit = artIter->second;
            hitVector.push_back(hit);

            const pandora::CartesianVector point(pCaloHit3D->GetPositionVector());
            double xyz[3] = { point.GetX(), point.GetY(), point.GetZ() };
            double dxdydz[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // TODO: Fill in errors
            double chi2(0.0);

            recob::SpacePoint newSpacePoint(xyz, dxdydz, chi2, spacePointCounter++);
            outputSpacePoints->push_back(newSpacePoint);

            util::CreateAssn(*this, evt, *(outputSpacePoints.get()), hitVector, *(outputSpacePointsToHits.get()));
            util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputSpacePoints.get()), *(outputParticlesToSpacePoints.get()),
                outputSpacePoints->size() - 1, outputSpacePoints->size());
        }

        // Build 2D Clusters   
        pandora::ClusterVector pandoraClusterVector(pPfo->GetClusterList().begin(), pPfo->GetClusterList().end());
        std::sort(pandoraClusterVector.begin(), pandoraClusterVector.end(), lar_content::LArClusterHelper::SortByNHits);
        HitVector pfoHits; // select all hits from this Pfo

        for (pandora::ClusterVector::const_iterator cIter = pandoraClusterVector.begin(), cIterEnd = pandoraClusterVector.end(); cIter != cIterEnd; ++cIter)
        {
            const pandora::Cluster *const pCluster = *cIter;

            if (pandora::TPC_3D == lar_content::LArClusterHelper::GetClusterHitType(pCluster))
                continue;

            pandora::CaloHitList pandoraHitList2D;
            pCluster->GetOrderedCaloHitList().GetCaloHitList(pandoraHitList2D);
            pandoraHitList2D.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
            pandora::CaloHitVector pandoraHitVector2D(pandoraHitList2D.begin(), pandoraHitList2D.end());
            std::sort(pandoraHitVector2D.begin(), pandoraHitVector2D.end(), lar_content::LArClusterHelper::SortByPosition);

            HitArray  hitArray;      // sort hits by drift volume
            HitList   isolatedHits;  // select isolated hits

            for (pandora::CaloHitVector::const_iterator hIter = pandoraHitVector2D.begin(), hIterEnd = pandoraHitVector2D.end(); hIter != hIterEnd; ++hIter)
            {
                const pandora::CaloHit *const pCaloHit = *hIter;

                const void *const pHitAddress(pCaloHit->GetParentCaloHitAddress());
                const intptr_t hitID_temp((intptr_t)(pHitAddress)); // TODO
                const int hitID((int)(hitID_temp));

                HitMap::const_iterator artIter = hitMap.find(hitID);

                if (hitMap.end() == artIter)
                    throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

                const art::Ptr<recob::Hit> hit = artIter->second;
                const geo::WireID wireID(hit->WireID());
                const unsigned int volID(100000 * wireID.Cryostat + wireID.TPC);

                hitArray[volID].push_back(hit);

                pfoHits.push_back(hit);
                if (pCaloHit->IsIsolated())
                    isolatedHits.insert(hit);
            }

            if (hitArray.empty())
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            for (HitArray::const_iterator hIter = hitArray.begin(), hIterEnd = hitArray.end(); hIter != hIterEnd; ++hIter)
            {
                const HitVector clusterHits(hIter->second);
                outputClusters->emplace_back(
                    LArPandoraHelper::BuildCluster(clusterCounter++, clusterHits, isolatedHits, ClusterParamAlgo)
                ); 

                util::CreateAssn(*this, evt, *(outputClusters.get()), clusterHits, *(outputClustersToHits.get()));
                util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputClusters.get()), *(outputParticlesToClusters.get()),
                    outputClusters->size() - 1, outputClusters->size());
            }
        }

        // Associate Vertex and Build Seeds (and Tracks)
        if (!pPfo->GetVertexList().empty())
        {
            if(pPfo->GetVertexList().size() != 1)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const pandora::Vertex *const pVertex = *(pPfo->GetVertexList().begin());

            ThreeDVertexMap::const_iterator vIter = vertexMap.find(pVertex);
            if (vertexMap.end() == vIter)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const unsigned int vtxElement(vIter->second);

            const pandora::CartesianVector vtxPos(pVertex->GetPosition());
            double pos[3] = { vtxPos.GetX(), vtxPos.GetY(), vtxPos.GetZ() };
            double posErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors

            util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputVertices.get()), *(outputParticlesToVertices.get()),
                vtxElement, vtxElement + 1);

            if (lar_content::LArPfoHelper::IsTrack(pPfo) && pPfo->GetMomentum().GetMagnitudeSquared() > std::numeric_limits<float>::epsilon())
            {
                const pandora::CartesianVector vtxDir(pPfo->GetMomentum().GetUnitVector()); 
                double dir[3]     = { vtxDir.GetX(), vtxDir.GetY(), vtxDir.GetZ() };
                double dirErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors

                recob::Seed newSeed(pos, dir, posErr, dirErr);
                outputSeeds->push_back(newSeed);

                util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputSeeds.get()), *(outputParticlesToSeeds.get()),
                    outputSeeds->size() - 1, outputSeeds->size());

                if (m_buildTracks)
                {
                    try
                    {
                        recob::Track newTrack(LArPandoraHelper::BuildTrack(trackCounter++, pPfo)); 
                        outputTracks->push_back(newTrack);

                        util::CreateAssn(*this, evt, *(outputTracks.get()), pfoHits, *(outputTracksToHits.get()));
                        util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputTracks.get()), *(outputParticlesToTracks.get()),
                            outputTracks->size() - 1, outputTracks->size());
                    }
                    catch (cet::exception &e)
                    {
                    }
                }
            }
        }

        //
        // TODO: Build Showers here
        //
    }

    mf::LogDebug("LArPandora") << "   Number of new particles: " << outputParticles->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new clusters: " << outputClusters->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new space points: " << outputSpacePoints->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new seeds: " << outputSeeds->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new vertices: " << outputVertices->size() << std::endl;

    if (m_buildTracks)
        mf::LogDebug("LArPandora") << "   Number of new tracks: " << outputTracks->size() << std::endl;

    if (m_buildShowers)
        mf::LogDebug("LArPandora") << "   Number of new showers: " << outputShowers->size() << std::endl;

    evt.put(std::move(outputParticles));
    evt.put(std::move(outputSpacePoints));
    evt.put(std::move(outputClusters));
    evt.put(std::move(outputSeeds));
    evt.put(std::move(outputVertices));

    evt.put(std::move(outputParticlesToSpacePoints));
    evt.put(std::move(outputParticlesToClusters));
    evt.put(std::move(outputParticlesToSeeds));
    evt.put(std::move(outputParticlesToVertices));
    evt.put(std::move(outputSpacePointsToHits));
    evt.put(std::move(outputClustersToHits));

    if (m_buildTracks)
    {
        evt.put(std::move(outputTracks));
        evt.put(std::move(outputParticlesToTracks));
        evt.put(std::move(outputTracksToHits));
    }

    if (m_buildShowers)
    {
        evt.put(std::move(outputShowers));
        evt.put(std::move(outputParticlesToShowers));
        evt.put(std::move(outputShowersToHits));
    }
}

} // namespace lar_pandora
