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

// Pandora includes
#include "Objects/ParticleFlowObject.h"

// Local includes (LArPandoraAlgorithms)
#include "LArContent.h"
#include "LArHelpers/LArPfoHelper.h"

// Local includes (LArPandoraInterface)
#include "LArPandoraParticleCreator.h"
#include "LArPandoraHelper.h"
#include "PFParticleSeed.h"

// System includes
#include <iostream>
#include <limits>

namespace lar_pandora {

LArPandoraParticleCreator::LArPandoraParticleCreator(fhicl::ParameterSet const &pset) : LArPandoraBase(pset)
{
    produces< std::vector<recob::PFParticle> >();
    produces< std::vector<recob::SpacePoint> >();
    produces< std::vector<recob::Cluster> >(); 
    produces< std::vector<recob::Seed> >();

    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Seed> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
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
    {
        mf::LogDebug("LArPandora") << "   No reconstructed particles for this event [return] " << std::endl;
        return;
    }


    // Set up ART outputs
    std::unique_ptr< std::vector<recob::PFParticle> > outputParticles( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::SpacePoint> > outputSpacePoints( new std::vector<recob::SpacePoint> );
    std::unique_ptr< std::vector<recob::Cluster> >    outputClusters( new std::vector<recob::Cluster> );
    std::unique_ptr< std::vector<recob::Seed> >       outputSeeds( new std::vector<recob::Seed> );

    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >    outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> >    outputParticlesToSeeds( new art::Assns<recob::PFParticle, recob::Seed> );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> >   outputSpacePointsToHits( new art::Assns<recob::SpacePoint, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> >      outputClustersToHits( new art::Assns<recob::Cluster, recob::Hit> );


    // Obtain a sorted vector of all output Pfos and their daughters
    pandora::PfoList connectedPfoList;
    lar_content::LArPfoHelper::GetAllConnectedPfos(concatenatedPfoList, connectedPfoList);
    pandora::PfoVector pfoVector(connectedPfoList.begin(), connectedPfoList.end());
    std::sort(pfoVector.begin(), pfoVector.end(), lar_content::LArPfoHelper::SortByNHits);

    typedef std::map< const pandora::ParticleFlowObject*, size_t> ThreeDParticleMap;
    typedef std::map< int, HitVector > HitArray;

    int spacePointCounter(0);
    int clusterCounter(0);
    size_t particleCounter(0);

    // Build map of Pandora particles
    ThreeDParticleMap particleMap;
    for (pandora::PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const pandora::ParticleFlowObject *const pPfo = *pIter;
        particleMap.insert( std::pair<const pandora::ParticleFlowObject*, size_t>(pPfo, particleCounter++) );
    }

    // Loop over Pandora particles
    for (pandora::PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const pandora::ParticleFlowObject *const pPfo = *pIter;
        const pandora::ClusterList &pfoClusterList = pPfo->GetClusterList();

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

        // Get Pandora 3D Hits
        pandora::CaloHitList pandoraHitList3D;
        for (pandora::ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const pandora::Cluster *const pCluster = *cIter;

            if (pandora::TPC_3D != lar_content::LArClusterHelper::GetClusterHitType(pCluster))
                continue;

            pCluster->GetOrderedCaloHitList().GetCaloHitList(pandoraHitList3D);
        }

        // Build Particle
        recob::PFParticle newParticle(pPfo->GetParticleId(), pfoIdCode, parentIdCode, daughterIdCodes);
        outputParticles->push_back(newParticle);

        // Build Space Points (TODO: Order 3D hits according to displacement along 3D cluster)
	pandora::CartesianPointList pointList;

        for (pandora::CaloHitList::const_iterator hIter = pandoraHitList3D.begin(), hIterEnd = pandoraHitList3D.end(); hIter != hIterEnd; ++hIter)
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
            pointList.push_back(point);

            double xyz[3] = { point.GetX(), point.GetY(), point.GetZ() };
            double dxdydz[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // TODO: Fill in errors
            double chi2(0.0);

            recob::SpacePoint newSpacePoint(xyz, dxdydz, chi2, spacePointCounter++);
            outputSpacePoints->push_back(newSpacePoint);

            util::CreateAssn(*this, evt, *(outputSpacePoints.get()), hitVector, *(outputSpacePointsToHits.get()));
            util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputSpacePoints.get()), *(outputParticlesToSpacePoints.get()),
                outputSpacePoints->size() - 1, outputSpacePoints->size());
        }

        // Build Seeds (TODO: Persist these into the data stream)
	if (lar_content::LArPfoHelper::IsTrack(pPfo) && !pointList.empty())
	{
	    const PFParticleSeed seed(pointList); // (TODO: Get this information from Pandora output)
            if (seed.IsInitialized())
	    {
	        for (int n=0; n<2; ++n)
		{
		    const pandora::CartesianVector vtxPos((0 == n) ? seed.GetInnerPosition() : seed.GetOuterPosition());
                    const pandora::CartesianVector vtxDir((0 == n) ? seed.GetInnerDirection() : seed.GetOuterDirection());
	        
	            double pos[3]     = { vtxPos.GetX(), vtxPos.GetY(), vtxPos.GetZ() };
                    double dir[3]     = { vtxDir.GetX(), vtxDir.GetY(), vtxDir.GetZ() };
                    double posErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors
                    double dirErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors
		
		    recob::Seed newSeed(pos, dir, posErr, dirErr);
                    outputSeeds->push_back(newSeed);

                    util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputSeeds.get()), *(outputParticlesToSeeds.get()),
                        outputSeeds->size() - 1, outputSeeds->size());
		}
	    }
	}

        // Build Clusters 
        for (pandora::ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const pandora::Cluster *const pCluster = *cIter;

            if (pandora::TPC_3D == lar_content::LArClusterHelper::GetClusterHitType(pCluster))
                continue;

            pandora::CaloHitList pandoraHitList2D;
            pCluster->GetOrderedCaloHitList().GetCaloHitList(pandoraHitList2D);
            pandoraHitList2D.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            HitArray hitArray;
            for (pandora::CaloHitList::const_iterator hIter = pandoraHitList2D.begin(), hIterEnd = pandoraHitList2D.end(); hIter != hIterEnd; ++hIter)
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
            }

            if (hitArray.empty())
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            for (HitArray::const_iterator hIter = hitArray.begin(), hIterEnd = hitArray.end(); hIter != hIterEnd; ++hIter)
	    {
	        const HitVector hitVector(hIter->second);
                recob::Cluster newCluster(LArPandoraHelper::BuildCluster(clusterCounter++, hitVector));
                outputClusters->push_back(newCluster);
	    
                util::CreateAssn(*this, evt, *(outputClusters.get()), hitVector, *(outputClustersToHits.get()));
                util::CreateAssn(*this, evt, *(outputParticles.get()), *(outputClusters.get()), *(outputParticlesToClusters.get()),
                    outputClusters->size() - 1, outputClusters->size());
	    }
        }
    }

    mf::LogDebug("LArPandora") << "   Number of new particles: " << outputParticles->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new clusters: " << outputClusters->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new space points: " << outputSpacePoints->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new seeds: " << outputSeeds->size() << std::endl;

    evt.put(std::move(outputParticles));
    evt.put(std::move(outputSpacePoints));
    evt.put(std::move(outputClusters));
    evt.put(std::move(outputSeeds));

    evt.put(std::move(outputParticlesToSpacePoints));
    evt.put(std::move(outputParticlesToClusters));
    evt.put(std::move(outputParticlesToSeeds));
    evt.put(std::move(outputSpacePointsToHits));
    evt.put(std::move(outputClustersToHits));
}

} // namespace lar_pandora
