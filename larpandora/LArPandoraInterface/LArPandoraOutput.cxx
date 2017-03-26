/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraOutput.cxx
 *
 *  @brief  Helper functions for processing outputs from pandora
 *
 */

#include "art/Framework/Core/EDProducer.h"
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::TPCID
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/PtrMaker.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "Api/PandoraApi.h"

#include "Objects/ParticleFlowObject.h"
#include "Objects/TrackState.h"
#include "Objects/Vertex.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <limits>

namespace lar_pandora
{

void LArPandoraOutput::ProduceArtOutput(const Settings &settings, const IdToHitMap &idToHitMap, art::Event &evt)
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ProduceArtOutput() *** " << std::endl;

    if (!settings.m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " Throwing exception - primary Pandora instance does not exist ";

    if (!settings.m_pProducer)
        throw cet::exception("LArPandora") << " Throwing exception - pointer to ART Producer module does not exist ";

    PandoraInstanceList pandoraInstanceList;
    const PandoraInstanceList &daughterInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(settings.m_pPrimaryPandora));

    if (settings.m_buildStitchedParticles || daughterInstances.empty())
    {
        pandoraInstanceList.push_back(settings.m_pPrimaryPandora);
    }
    else
    {
        pandoraInstanceList.insert(pandoraInstanceList.end(), daughterInstances.begin(), daughterInstances.end());
    }

    pandora::PfoList concatenatedPfoList;
    for (const pandora::Pandora *const pPandora : pandoraInstanceList)
    {
        const pandora::PfoList *pPfoList(nullptr);
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pPandora, pPfoList));
        concatenatedPfoList.insert(concatenatedPfoList.end(), pPfoList->begin(), pPfoList->end());
    }

    if (concatenatedPfoList.empty())
        mf::LogDebug("LArPandora") << "   Warning: No reconstructed particles for this event " << std::endl;

    // Set up ART outputs from RecoBase and AnalysisBase
    std::unique_ptr< std::vector<recob::PFParticle> > outputParticles( new std::vector<recob::PFParticle> );
    std::unique_ptr< std::vector<recob::SpacePoint> > outputSpacePoints( new std::vector<recob::SpacePoint> );
    std::unique_ptr< std::vector<recob::Cluster> >    outputClusters( new std::vector<recob::Cluster> );
    std::unique_ptr< std::vector<recob::Seed> >       outputSeeds( new std::vector<recob::Seed> );
    std::unique_ptr< std::vector<recob::Vertex> >     outputVertices( new std::vector<recob::Vertex> );
    std::unique_ptr< std::vector<recob::Track> >      outputTracks( new std::vector<recob::Track> );
    std::unique_ptr< std::vector<recob::Shower> >     outputShowers( new std::vector<recob::Shower> );
    std::unique_ptr< std::vector<recob::PCAxis> >     outputPCAxes( new std::vector<recob::PCAxis> );
    std::unique_ptr< std::vector<anab::T0> >          outputT0s( new std::vector<anab::T0> );

    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >    outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed> >       outputParticlesToSeeds( new art::Assns<recob::PFParticle, recob::Seed> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> >     outputParticlesToVertices( new art::Assns<recob::PFParticle, recob::Vertex> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Track> >      outputParticlesToTracks( new art::Assns<recob::PFParticle, recob::Track> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower> >     outputParticlesToShowers( new art::Assns<recob::PFParticle, recob::Shower> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> >     outputParticlesToPCAxes( new art::Assns<recob::PFParticle, recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::Track, recob::Hit> >             outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit> >            outputShowersToHits( new art::Assns<recob::Shower, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Shower, recob::PCAxis> >         outputShowersToPCAxes( new art::Assns<recob::Shower, recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> >        outputSpacePointsToHits( new art::Assns<recob::SpacePoint, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> >           outputClustersToHits( new art::Assns<recob::Cluster, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Seed, recob::Hit> >              outputSeedsToHits( new art::Assns<recob::Seed, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Track, anab::T0> >               outputTracksToT0s( new art::Assns<recob::Track, anab::T0> );

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    cluster::StandardClusterParamsAlg ClusterParamAlgo;

    // Obtain a sorted vector of all output Pfos and their daughters
    pandora::PfoList connectedPfoList;
    lar_content::LArPfoHelper::GetAllConnectedPfos(concatenatedPfoList, connectedPfoList);

    pandora::PfoVector pfoVector(connectedPfoList.begin(), connectedPfoList.end());
    std::sort(pfoVector.begin(), pfoVector.end(), lar_content::LArPfoHelper::SortByNHits);

    int vertexCounter(0);
    int spacePointCounter(0);
    int clusterCounter(0);
    size_t clusterHitAssnCounter(0);
    int trackCounter(0);
    size_t particleCounter(0);
    int t0Counter(0);

    // Build maps of Pandora particles and Pandora vertices
    pandora::VertexVector vertexVector;
    ThreeDParticleMap particleMap;
    ThreeDVertexMap vertexMap;

    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
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

    auto const& geom = lar::providerFrom<geo::Geometry>();

    // Loop over Pandora vertices and build recob::Vertices
    for (const pandora::Vertex *const pVertex : vertexVector)
    {
        ThreeDVertexMap::const_iterator iter = vertexMap.find(pVertex);
        if (vertexMap.end() == iter)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const pandora::CartesianVector vtxPos(pVertex->GetPosition());
        double pos[3] = { vtxPos.GetX(), vtxPos.GetY(), vtxPos.GetZ() };

        recob::Vertex newVertex(pos, vertexCounter++);
        outputVertices->push_back(newVertex);
    }

    lar::PtrMaker<recob::Shower> makeShowerPtr(evt, *(settings.m_pProducer));
    lar::PtrMaker<recob::PCAxis> makePCAxisPtr(evt, *(settings.m_pProducer));
    lar::PtrMaker<recob::PFParticle> makePfoPtr(evt, *(settings.m_pProducer));
    lar::PtrMaker<recob::Cluster> makeClusterPtr(evt, *(settings.m_pProducer));

    // Loop over Pandora particles and build recob::PFParticles
    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
        // Get Pfo ID
        ThreeDParticleMap::const_iterator iter = particleMap.find(pPfo);
        if (particleMap.end() == iter)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

        const size_t pfoIdCode(iter->second);

        // Get Pfo Parents
        size_t parentIdCode(recob::PFParticle::kPFParticlePrimary);
        const pandora::PfoList &parentList(pPfo->GetParentPfoList());

        if (!parentList.empty())
        {
            if (parentList.size() != 1)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            ThreeDParticleMap::const_iterator parentIdIter = particleMap.find(*parentList.begin());
            if (particleMap.end() == parentIdIter)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            parentIdCode = parentIdIter->second;
        }

        // Get Pfo Daughters
        std::vector<size_t> daughterIdCodes;
        pandora::PfoVector daughterPfoVector(pPfo->GetDaughterPfoList().begin(), pPfo->GetDaughterPfoList().end());
        std::sort(daughterPfoVector.begin(), daughterPfoVector.end(), lar_content::LArPfoHelper::SortByNHits);

        for (const pandora::ParticleFlowObject *const pDaughterPfo : daughterPfoVector)
        {
            ThreeDParticleMap::const_iterator daughterIdIter = particleMap.find(pDaughterPfo);
            if (particleMap.end() == daughterIdIter)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const size_t daughterIdCode(daughterIdIter->second);
            daughterIdCodes.push_back(daughterIdCode);
        }

        // Build Particle
        recob::PFParticle newParticle(pPfo->GetParticleId(), pfoIdCode, parentIdCode, daughterIdCodes);
        outputParticles->push_back(newParticle);

        // Store 2D Hits (from either clusters or space points)
        HitVector particleHitsFromClusters, particleHitsFromSpacePoints;

        // Build 2D Clusters
        pandora::ClusterVector pandoraClusterVector(pPfo->GetClusterList().begin(), pPfo->GetClusterList().end());
        std::sort(pandoraClusterVector.begin(), pandoraClusterVector.end(), lar_content::LArClusterHelper::SortByNHits);

        int iClusterCounter = clusterCounter;                   // Record placeholders for clusters and hits in preparation for shower-building
        size_t iClusterHitAssnCounter = clusterHitAssnCounter;  //

        for (const pandora::Cluster *const pCluster : pandoraClusterVector)
        {
            if (pandora::TPC_3D == lar_content::LArClusterHelper::GetClusterHitType(pCluster))
                continue;

            pandora::CaloHitList pandoraHitList2D;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraHitList2D);
            pandoraHitList2D.insert(pandoraHitList2D.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            pandora::CaloHitVector pandoraHitVector2D(pandoraHitList2D.begin(), pandoraHitList2D.end());
            std::sort(pandoraHitVector2D.begin(), pandoraHitVector2D.end(), lar_content::LArClusterHelper::SortHitsByPosition);

            HitArray  hitArray;      // sort hits by drift volume
            HitList   isolatedHits;  // select isolated hits

            for (const pandora::CaloHit *const pCaloHit2D : pandoraHitVector2D)
            {
                const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, pCaloHit2D);
                particleHitsFromClusters.push_back(hit);

                const geo::WireID wireID(hit->WireID());
                const unsigned int volID(100000 * wireID.Cryostat + wireID.TPC);

                hitArray[volID].push_back(hit);

                if (pCaloHit2D->IsIsolated())
                    isolatedHits.insert(hit);
            }

            if (hitArray.empty())
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            for (const HitArray::value_type &hitArrayEntry : hitArray)
            {
                const HitVector &clusterHits(hitArrayEntry.second);
                outputClusters->emplace_back(LArPandoraOutput::BuildCluster(clusterCounter++, clusterHits, isolatedHits, ClusterParamAlgo));
                clusterHitAssnCounter += clusterHits.size();

                util::CreateAssn(*(settings.m_pProducer), evt, *(outputClusters.get()), clusterHits, *(outputClustersToHits.get()));
                util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputClusters.get()), *(outputParticlesToClusters.get()),
                    outputClusters->size() - 1, outputClusters->size());

                LOG_DEBUG("LArPandora") << "Stored cluster ID="
                  << outputClusters->back().ID()
                  << " (#" << (outputClusters->size() - 1)
                  << ") with " << clusterHits.size() << " hits";
            }
        }

        // Build 3D SpacePoints and calculate T0 from shift in 2D hits
        pandora::CaloHitList pandoraHitList3D;
        lar_content::LArPfoHelper::GetCaloHits(pPfo, pandora::TPC_3D, pandoraHitList3D);

        pandora::CaloHitVector pandoraHitVector3D(pandoraHitList3D.begin(), pandoraHitList3D.end());
        std::sort(pandoraHitVector3D.begin(), pandoraHitVector3D.end(), lar_content::LArClusterHelper::SortHitsByPosition);

        double sumT(0.0), sumN(0.0);

        for (const pandora::CaloHit *const pCaloHit3D : pandoraHitVector3D)
        {
            if (pandora::TPC_3D != pCaloHit3D->GetHitType())
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const pandora::CaloHit *const pCaloHit2D = static_cast<const pandora::CaloHit*>(pCaloHit3D->GetParentAddress());

            const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, pCaloHit2D);
            particleHitsFromSpacePoints.push_back(hit);

            HitVector spacePointHits;
            spacePointHits.push_back(hit);

            // ATTN: We assume that the 2D Pandora hits have been shifted
            sumT += LArPandoraOutput::CalculateT0(hit, pCaloHit2D);
            sumN += 1.0;

            outputSpacePoints->emplace_back(LArPandoraOutput::BuildSpacePoint(spacePointCounter++, pCaloHit3D));

            util::CreateAssn(*(settings.m_pProducer), evt, *(outputSpacePoints.get()), spacePointHits, *(outputSpacePointsToHits.get()));
            util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputSpacePoints.get()), *(outputParticlesToSpacePoints.get()),
                outputSpacePoints->size() - 1, outputSpacePoints->size());
        }

        const double T0((sumN > 0.0 && std::fabs(sumT) > sumN) ? (sumT / sumN) : 0.0);

        // Associate Vertex and Build High-Level Objects
        if (!pPfo->GetVertexList().empty())
        {
            if(pPfo->GetVertexList().size() != 1)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const pandora::Vertex *const pVertex = *(pPfo->GetVertexList().begin());

            ThreeDVertexMap::const_iterator iter = vertexMap.find(pVertex);
            if (vertexMap.end() == iter)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

            const unsigned int vtxElement(iter->second);
            util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputVertices.get()), *(outputParticlesToVertices.get()),
                vtxElement, vtxElement + 1);

            // Build Seeds, Tracks, T0s
            if (lar_content::LArPfoHelper::IsTrack(pPfo) && pPfo->GetMomentum().GetMagnitudeSquared() > std::numeric_limits<float>::epsilon())
            {
                const lar_content::LArTrackPfo *const pLArTrackPfo = dynamic_cast<const lar_content::LArTrackPfo*>(pPfo);

                if (!pLArTrackPfo)
                {
                    mf::LogDebug("LArPandoraOutput") << " LArPandoraOutput::BuildTrack --- input pfo is track-like but is not a LArTrackPfo ";
                    continue;
                }

                const lar_content::LArTrackStateVector &trackStateVector = pLArTrackPfo->m_trackStateVector;

                if (trackStateVector.size() < settings.m_minTrajectoryPoints)
                {
                    mf::LogDebug("LArPandoraOutput") << " LArPandoraOutput::BuildTrack --- Insufficient input trajectory points to build track ";
                    continue;
                }

                for (const lar_content::LArTrackState &nextPoint : trackStateVector)
                {
                    const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, nextPoint.GetCaloHit());

                    HitVector seedHits;
                    seedHits.push_back(hit);

                    outputSeeds->emplace_back(LArPandoraOutput::BuildSeed(nextPoint));

                    util::CreateAssn(*(settings.m_pProducer), evt, *(outputSeeds.get()), seedHits, *(outputSeedsToHits.get()));
                    util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputSeeds.get()), *(outputParticlesToSeeds.get()),
                        outputSeeds->size() - 1, outputSeeds->size());
                }

                if (!settings.m_buildTracks)
                    continue;

                // Building track objects
                outputTracks->emplace_back(LArPandoraOutput::BuildTrack(trackCounter++, &trackStateVector, idToHitMap));

                util::CreateAssn(*(settings.m_pProducer), evt, *(outputTracks.get()), particleHitsFromSpacePoints, *(outputTracksToHits.get()));
                util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputTracks.get()), *(outputParticlesToTracks.get()), outputTracks->size() - 1, outputTracks->size());

                // Output T0 objects [arguments are:  time (nanoseconds);  trigger type (3 for TPC stitching!);  track ID code;  T0 ID code]
                // ATTN: T0 values are currently calculated in nanoseconds relative to the trigger offset. Only non-zero values are outputted.
                if (std::fabs(T0) > 0.0)
                {
                    outputT0s->emplace_back(anab::T0(T0, 3, outputTracks->back().ID(), t0Counter++));
                    util::CreateAssn(*(settings.m_pProducer), evt, *(outputTracks.get()), *(outputT0s.get()), *(outputTracksToT0s.get()), outputT0s->size() - 1, outputT0s->size());
                }
            }

            // Build Showers, PCAxes
            else if (lar_content::LArPfoHelper::IsShower(pPfo))
            {
                const lar_content::LArShowerPfo *const pLArShowerPfo = dynamic_cast<const lar_content::LArShowerPfo*>(pPfo);

                if (!pLArShowerPfo)
                {
                    mf::LogDebug("LArPandoraOutput") << " LArPandoraOutput::BuildShower --- input pfo is shower-like but is not a LArShowerPfo ";
                    continue;
                }

                if (!settings.m_buildShowers)
                    continue;

                // TODO - If possible, we should try to move some of the shower-building code below into the BuildShower method

                // if this assertion fails, we have a shower with no associated clusters:
                assert((size_t) iClusterCounter < outputClusters->size());
                assert(iClusterHitAssnCounter < outputClustersToHits->size());

                std::vector<double> showerE; // will be empty if no energy algorithm was requested
                if(settings.m_showerEnergyAlg) {

                    // we expect all the clusters to be within this TPC ... ATTN: Will this work for DUNE?
                    geo::TPCID refTPC = (*outputClusters)[iClusterCounter].Plane();

                    // prepare the energies vector with one entry per plane
                    // (we get the total number of planes of the TPC  the cluster is in from geometry)
                    // and initialize them to a ridiculously negative number to start with
                    showerE.resize(
                        geom->TPC(refTPC).Nplanes(),
                        std::numeric_limits<double>::lowest()
                        );

                    size_t const nClusters = outputClusters->size() - iClusterCounter;

                    LOG_DEBUG("LArPandora")
                        << nClusters << " clusters for shower #" << outputShowers->size();
                    if ( nClusters > showerE.size() ) {
                        // not fun, but we push through
                        mf::LogError("LArPandora") << nClusters << " clusters for "
                            << showerE.size() << " wire planes!";
                    }

                    // go through the new clusters
                    // iterator to the first association of an hit to a new cluster:
                    auto beginHit = outputClustersToHits->begin() + iClusterHitAssnCounter;
                    for ( size_t iCluster = iClusterCounter; iCluster < outputClusters->size(); ++iCluster ) {
                        auto const& cluster = (*outputClusters)[iCluster];

                        auto const clusterPlaneID = cluster.Plane();
                        if (refTPC != clusterPlaneID) {
                            throw cet::exception("LArPandora")
                              << "Clusters for shower #" << outputShowers->size()
                              << " are expected on TPC " << std::string(refTPC)
                              << " but cluster ID=" << cluster.ID()
                              << " is on plane " + std::string(clusterPlaneID)
                              ;
                        }

                        //
                        // collect back the hits associated to this cluster (via art pointer)
                        //
                        std::vector<art::Ptr<recob::Hit>> clusterHits;
                        auto const clusterPtr = makeClusterPtr(iCluster);
                        auto endHit = beginHit;
                        while (endHit != outputClustersToHits->end()) {
                            if (endHit->first != clusterPtr) break;
                            clusterHits.push_back(endHit->second);
                            ++endHit;
                        } // while
                        LOG_TRACE("LArPandora")
                            << "  " << clusterHits.size() << " hits for cluster ID="
                            << cluster.ID() << " (#" << iCluster << ")";
                        if (clusterHits.empty()) {
                            // this is likely an error in the logic of this algorithm
                            throw cet::exception("LArPandora")
                                << "LArPandoraOutput::ProduceArtOutput(): no hits associated with a cluster!?";
                        }

                        //
                        // compute the energy for this cluster
                        //
                        double const E = settings.m_showerEnergyAlg->CalculateClusterEnergy
                            (cluster, clusterHits);
                        beginHit = endHit; // next iteration, we start from here

                        //
                        // store the energy in the cell pertaining the cluster plane
                        //
                        auto const planeNo = cluster.Plane().Plane;
                        if (showerE[planeNo] >= 0.) {
                            LOG_WARNING("LArPandora")
                                << "Warning! two or more clusters share plane "
                                << cluster.Plane() << "! (the last with energy " << E
                                << ", the previous " << showerE[planeNo] << " GeV)";
                        }
                        showerE[planeNo] = E;
                        LOG_TRACE("LArPandora") << "  cluster energy: " << E
                          << " GeV (plane: " << cluster.Plane() << ")";

                    } // for new clusters
                } // if shower energy

                // Save showers
                outputShowers->emplace_back(LArPandoraOutput::BuildShower(pLArShowerPfo, showerE));
                outputPCAxes->emplace_back(LArPandoraOutput::BuildShowerPCA(pLArShowerPfo));
                outputShowers->back().set_id(outputShowers->size()); // 1-based sequence

                outputParticlesToShowers->addSingle(
                    makePfoPtr(outputParticles->size() - 1),   // index of the PFO we just made
                    makeShowerPtr(outputShowers->size() - 1)   // index of the shower we just made
                    );
                outputParticlesToPCAxes->addSingle(makePfoPtr(outputParticles->size() - 1), makePCAxisPtr(outputPCAxes->size() - 1));
                outputShowersToPCAxes->addSingle(makeShowerPtr(outputShowers->size() - 1), makePCAxisPtr(outputPCAxes->size() - 1));

                // Save associations between showers and hits
                util::CreateAssn(*(settings.m_pProducer), evt, *(outputShowers.get()), particleHitsFromSpacePoints, *(outputShowersToHits.get()));
            }
        }
    } // for each reconstructed particle flow

    mf::LogDebug("LArPandora") << "   Number of new particles: " << outputParticles->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new clusters: " << outputClusters->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new space points: " << outputSpacePoints->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new seeds: " << outputSeeds->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new vertices: " << outputVertices->size() << std::endl;

    if (settings.m_buildTracks)
    {
        mf::LogDebug("LArPandora") << "   Number of new tracks: " << outputTracks->size() << std::endl;
        mf::LogDebug("LArPandora") << "   Number of new T0s: " << outputT0s->size() << std::endl;
    }

    if (settings.m_buildShowers)
    {
        mf::LogDebug("LArPandora") << "   Number of new showers: " << outputShowers->size() << std::endl;
        mf::LogDebug("LArPandora") << "   Number of new pcaxes:  " << outputPCAxes->size() << std::endl;
    }

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
    evt.put(std::move(outputSeedsToHits));

    if (settings.m_buildTracks)
    {
        evt.put(std::move(outputTracks));
        evt.put(std::move(outputParticlesToTracks));
        evt.put(std::move(outputTracksToHits));
        evt.put(std::move(outputT0s));
        evt.put(std::move(outputTracksToT0s));
    }

    if (settings.m_buildShowers)
    {
        evt.put(std::move(outputShowers));
        evt.put(std::move(outputParticlesToShowers));
        evt.put(std::move(outputShowersToHits));
        evt.put(std::move(outputPCAxes));
        evt.put(std::move(outputParticlesToPCAxes));
        evt.put(std::move(outputShowersToPCAxes));
    }

    mf::LogDebug("LArPandora") << " *** LArPandora::ProduceArtOutput() [DONE!] *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Cluster LArPandoraOutput::BuildCluster(const int id, const HitVector &hitVector, const HitList &isolatedHits, cluster::ClusterParamsAlgBase &algo)
{
    mf::LogDebug("LArPandora") << "   Building Cluster [" << id << "], Number of hits = " << hitVector.size() << std::endl;

    if (hitVector.empty())
        throw cet::exception("LArPandora") << " LArPandoraOutput::BuildCluster --- No input hits were provided ";

    // Fill list of cluster properties
    geo::View_t view(geo::kUnknown);
    geo::PlaneID planeID;

    double startWire(+std::numeric_limits<float>::max()), sigmaStartWire(0.0);
    double startTime(+std::numeric_limits<float>::max()), sigmaStartTime(0.0);
    double endWire(-std::numeric_limits<float>::max()), sigmaEndWire(0.0);
    double endTime(-std::numeric_limits<float>::max()), sigmaEndTime(0.0);

    std::vector<recob::Hit const*> hits_for_params;
    hits_for_params.reserve(hitVector.size());

    for (const art::Ptr<recob::Hit> &hit : hitVector)
    {
        const double thisWire(hit->WireID().Wire);
        const double thisWireSigma(0.5);
        const double thisTime(hit->PeakTime());
        const double thisTimeSigma(double(2.*hit->RMS()));
        const geo::View_t thisView(hit->View());
        const geo::PlaneID thisPlaneID(hit->WireID().planeID());

        if (geo::kUnknown == view)
        {
            view = thisView;
            planeID = thisPlaneID;
        }

        if (!(thisView == view && thisPlaneID == planeID))
        {
            throw cet::exception("LArPandora") << " LArPandoraOutput::BuildCluster --- Input hits have inconsistent plane IDs ";
        }

        hits_for_params.push_back(&*hit);

        if (isolatedHits.count(hit))
            continue;

        if (thisWire < startWire || (thisWire == startWire && thisTime < startTime))
        {
            startWire = thisWire;
            sigmaStartWire = thisWireSigma;
            startTime = thisTime;
            sigmaStartTime = thisTimeSigma;
        }

        if (thisWire > endWire || (thisWire == endWire && thisTime > endTime))
        {
            endWire = thisWire;
            sigmaEndWire = thisWireSigma;
            endTime = thisTime;
            sigmaEndTime = thisTimeSigma;
        }

    }

    // feed the algorithm with all the cluster hits
    algo.SetHits(hits_for_params);

    // create the recob::Cluster directly in the vector
    return cluster::ClusterCreator(
      algo,                  // algo
      startWire,             // start_wire
      sigmaStartWire,        // sigma_start_wire
      startTime,             // start_tick
      sigmaStartTime,        // sigma_start_tick
      endWire,               // end_wire
      sigmaEndWire,          // sigma_end_wire
      endTime,               // end_tick
      sigmaEndTime,          // sigma_end_tick
      id,                    // ID
      view,                  // view
      planeID,               // plane
      recob::Cluster::Sentry // sentry
      ).move();
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Seed LArPandoraOutput::BuildSeed(const lar_content::LArTrackState &trackState)
{
    double pos[3]     = { trackState.GetPosition().GetX(), trackState.GetPosition().GetY(), trackState.GetPosition().GetZ() };
    double dir[3]     = { trackState.GetDirection().GetX(), trackState.GetDirection().GetY(), trackState.GetDirection().GetZ() };
    double posErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors
    double dirErr[3]  = { 0.0, 0.0, 0.0 };  // TODO: Fill in errors

    return recob::Seed(pos, dir, posErr, dirErr);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Track LArPandoraOutput::BuildTrack(const int id, const lar_content::LArTrackStateVector *const pTrackStateVector, const IdToHitMap &idToHitMap)
{
    mf::LogDebug("LArPandora") << "   Building Track [" << id << "], Number of trajectory points = " << pTrackStateVector->size() << std::endl;

    if (pTrackStateVector->empty())
        throw cet::exception("LArPandora") << " LArPandoraOutput::BuildTrack --- No input trajectory points were provided ";

    // Fill list of track properties
    std::vector<TVector3>               xyz;
    std::vector<TVector3>               pxpypz;
    std::vector< std::vector<double> >  dQdx(3);
    std::vector<double>                 momentum = std::vector<double>(2, util::kBogusD);

    // Loop over trajectory points
    for (const lar_content::LArTrackState &nextPoint : *pTrackStateVector)
    {
        const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, nextPoint.GetCaloHit());
        const geo::View_t hit_View(hit->View());

        const TVector3 trackPosition(nextPoint.GetPosition().GetX(), nextPoint.GetPosition().GetY(), nextPoint.GetPosition().GetZ());
        const TVector3 trackDirection(nextPoint.GetDirection().GetX(), nextPoint.GetDirection().GetY(), nextPoint.GetDirection().GetZ());
        const double trackdQdx(LArPandoraOutput::CalculatedQdL(hit, trackPosition, trackDirection));

        const double dQdxU((geo::kU == hit_View) ? trackdQdx : 0.0);
        const double dQdxV((geo::kV == hit_View) ? trackdQdx : 0.0);
        const double dQdxW((geo::kW == hit_View) ? trackdQdx : 0.0);

        xyz.push_back(trackPosition);
        pxpypz.push_back(trackDirection);

        dQdx.at(geo::kU).push_back(dQdxU); dQdx.at(geo::kV).push_back(dQdxV); dQdx.at(geo::kW).push_back(dQdxW);
    }

    // Return a new recob::Track object (of the Bezier variety)
    return recob::Track(xyz, pxpypz, dQdx, momentum, id);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Shower LArPandoraOutput::BuildShower(const lar_content::LArShowerPfo *const pLArShowerPfo, const std::vector<double>& totalEnergy)
{
    const pandora::CartesianVector &showerLength(pLArShowerPfo->GetShowerLength());
    const pandora::CartesianVector &showerDirection(pLArShowerPfo->GetShowerDirection());
    const pandora::CartesianVector &showerVertex(pLArShowerPfo->GetShowerVertex());

    const float length(showerLength.GetX());
    const float openingAngle(pLArShowerPfo->GetShowerOpeningAngle());
    const TVector3 direction(showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ());
    const TVector3 vertex(showerVertex.GetX(), showerVertex.GetY(), showerVertex.GetZ());

    // TODO
    const TVector3 directionErr;
    const TVector3 vertexErr;
    const std::vector<double> totalEnergyErr;
    const std::vector<double> dEdx;
    const std::vector<double> dEdxErr;
    const int bestplane(0);

    return recob::Shower(direction, directionErr, vertex, vertexErr, totalEnergy, totalEnergyErr, dEdx, dEdxErr, bestplane, util::kBogusI, length, openingAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::PCAxis LArPandoraOutput::BuildShowerPCA(const lar_content::LArShowerPfo *const pLArShowerPfo)
{
    const pandora::CartesianVector &showerCentroid(pLArShowerPfo->GetShowerCentroid());
    const pandora::CartesianVector &showerDirection(pLArShowerPfo->GetShowerDirection());
    const pandora::CartesianVector &showerSecondaryVector(pLArShowerPfo->GetShowerSecondaryVector());
    const pandora::CartesianVector &showerTertiaryVector(pLArShowerPfo->GetShowerTertiaryVector());
    const pandora::CartesianVector &showerEigenValues(pLArShowerPfo->GetShowerEigenValues());

    const bool svdOK(true); ///< SVD Decomposition was successful
    const double eigenValues[3] = {showerEigenValues.GetX(), showerEigenValues.GetY(), showerEigenValues.GetZ()}; ///< Eigen values from SVD decomposition
    const double avePosition[3] = {showerCentroid.GetX(), showerCentroid.GetY(), showerCentroid.GetZ()}; ///< Average position of hits fed to PCA

    std::vector< std::vector<double> > eigenVecs = { /// The three principle axes
        { showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ() },
        { showerSecondaryVector.GetX(), showerSecondaryVector.GetY(), showerSecondaryVector.GetZ() },
        { showerTertiaryVector.GetX(), showerTertiaryVector.GetY(), showerTertiaryVector.GetZ() }
    };

    // TODO
    const int numHitsUsed(100); ///< Number of hits in the decomposition, not yet ready
    const double aveHitDoca(0.); ///< Average doca of hits used in PCA, not ready yet
    const size_t iD(util::kBogusI); ///< Axis ID, not ready yet

    return recob::PCAxis(svdOK, numHitsUsed, eigenValues, eigenVecs, avePosition, aveHitDoca, iD);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::SpacePoint LArPandoraOutput::BuildSpacePoint(const int id, const pandora::CaloHit *const pCaloHit)
{
    if (pandora::TPC_3D != pCaloHit->GetHitType())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const pandora::CartesianVector point(pCaloHit->GetPositionVector());
    double xyz[3] = { point.GetX(), point.GetY(), point.GetZ() };
    double dxdydz[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // TODO: Fill in the error matrix
    double chi2(0.0);

    return recob::SpacePoint(xyz, dxdydz, chi2, id);
}

//------------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Hit> LArPandoraOutput::GetHit(const IdToHitMap &idToHitMap, const pandora::CaloHit *const pCaloHit)
{
    const void *const pHitAddress(pCaloHit->GetParentAddress());
    const intptr_t hitID_temp((intptr_t)(pHitAddress));
    const int hitID((int)(hitID_temp));

    IdToHitMap::const_iterator artIter = idToHitMap.find(hitID);

    if (idToHitMap.end() == artIter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    return artIter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArPandoraOutput::CalculateT0(const art::Ptr<recob::Hit> hit, const pandora::CaloHit *const pCaloHit)
{
    art::ServiceHandle<geo::Geometry> theGeometry;
    auto const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();

    const geo::WireID hit_WireID(hit->WireID());
    const geo::TPCGeo &theTpc = theGeometry->Cryostat(hit_WireID.Cryostat).TPC(hit_WireID.TPC);

    // Calculate shift in x position between input and output hits
    const double input_xpos_cm(theDetector->ConvertTicksToX(hit->PeakTime(), hit_WireID.Plane, hit_WireID.TPC, hit_WireID.Cryostat));
    const double output_xpos_dm(pCaloHit->GetPositionVector().GetX());
    const double x0_cm(output_xpos_dm - input_xpos_cm);

    // The ingredients for the T0 calculation all come from the detector properties service
    const double dir((theTpc.DriftDirection() == geo::kNegX) ? 1.0 : -1.0);
    const double cm_per_tick(theDetector->GetXTicksCoefficient());
    const double ns_per_tick(theDetector->SamplingRate());

    // This calculation should give the T0 in nanoseconds relative to the initial 2D hit
    return (- dir * x0_cm * ns_per_tick / cm_per_tick);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArPandoraOutput::CalculatedQdL(const art::Ptr<recob::Hit> hit, const TVector3&, const TVector3 &trackDirection)
{
    // Extract wire pitch and wire direction from geometry
    art::ServiceHandle<geo::Geometry> theGeometry;
    const geo::WireID wireID(hit->WireID());
    const geo::WireGeo &wireGeo = theGeometry->Cryostat(wireID.Cryostat).TPC(wireID.TPC).Plane(wireID.Plane).Wire(wireID.Wire);
    const TVector3 driftDirection(1.0, 0.0, 0.0);
    const TVector3 wireDirection(wireGeo.Direction());
    const TVector3 wireAxis(wireDirection.Cross(driftDirection));
    const double wirePitch(theGeometry->WirePitch(hit->View()));

    // Calculate dQ/dL by projecting track direction onto wire
    const float cosTheta(std::fabs(wireAxis.Dot(trackDirection)));
    const float inverse_dL(cosTheta / wirePitch);
    const float dL((inverse_dL > std::numeric_limits<float>::epsilon()) ? (1.0 / inverse_dL) : std::numeric_limits<float>::max());
    const float dQ(hit->Integral());

    return (dQ/dL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraOutput::Settings::Settings() :
    m_pPrimaryPandora(nullptr),
    m_pProducer(nullptr),
    m_buildTracks(true),
    m_minTrajectoryPoints(2),
    m_buildShowers(true),
    m_buildStitchedParticles(false),
    m_showerEnergyAlg(nullptr)
{
}

} // namespace lar_pandora
