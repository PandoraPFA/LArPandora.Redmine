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
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "Api/PandoraApi.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
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
        throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- primary Pandora instance does not exist ";

    if (!settings.m_pProducer)
        throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- pointer to ART Producer module does not exist ";

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
    std::unique_ptr< std::vector<recob::Vertex> >     outputVertices( new std::vector<recob::Vertex> );
    std::unique_ptr< std::vector<anab::T0> >          outputT0s( new std::vector<anab::T0> );

    std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint> > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> >    outputParticlesToClusters( new art::Assns<recob::PFParticle, recob::Cluster> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> >     outputParticlesToVertices( new art::Assns<recob::PFParticle, recob::Vertex> );
    std::unique_ptr< art::Assns<recob::PFParticle, anab::T0> >          outputParticlesToT0s( new art::Assns<recob::PFParticle, anab::T0> );
    std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit> >        outputSpacePointsToHits( new art::Assns<recob::SpacePoint, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> >           outputClustersToHits( new art::Assns<recob::Cluster, recob::Hit> );

    // prepare the algorithm to compute the cluster characteristics;
    // we use the "standard" one here; configuration would happen here,
    // but we are using the default configuration for that algorithm
    cluster::StandardClusterParamsAlg ClusterParamAlgo;

    size_t particleCounter(0), vertexCounter(0), spacePointCounter(0), clusterCounter(0), t0Counter(0);

    // Obtain a sorted vector of all output Pfos and their daughters
    pandora::PfoList connectedPfoList;
    lar_content::LArPfoHelper::GetAllConnectedPfos(concatenatedPfoList, connectedPfoList);

    pandora::PfoVector pfoVector(connectedPfoList.begin(), connectedPfoList.end());
    std::sort(pfoVector.begin(), pfoVector.end(), lar_content::LArPfoHelper::SortByNHits);

    // Build maps of pandora::Pfos and build recob::vertices
    ThreeDParticleMap particleMap;
    ThreeDVertexMap vertexMap;

    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
        particleMap.insert( std::pair<const pandora::ParticleFlowObject*, size_t>(pPfo, particleCounter++) );

        if (!pPfo->GetVertexList().empty())
        {
            if(pPfo->GetVertexList().size() != 1)
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- this particle has multiple interaction vertices ";

            const pandora::Vertex *const pVertex(pPfo->GetVertexList().front());

            if (vertexMap.end() != vertexMap.find(pVertex))
                continue;

            double pos[3] = {pVertex->GetPosition().GetX(), pVertex->GetPosition().GetY(), pVertex->GetPosition().GetZ()};
            outputVertices->emplace_back(recob::Vertex(pos, vertexCounter++));
            vertexMap.insert(std::pair<const pandora::Vertex*, unsigned int>(pVertex, vertexCounter - 1));
        }
    }

    // Loop over pandora::Pfos and build recob::PFParticles
    for (const pandora::ParticleFlowObject *const pPfo : pfoVector)
    {
        // Get Pfo ID
        ThreeDParticleMap::const_iterator iter = particleMap.find(pPfo);
        if (particleMap.end() == iter)
            throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- found an unassociated particle";

        const size_t pfoIdCode(iter->second);

        // Get Pfo Parents
        size_t parentIdCode(recob::PFParticle::kPFParticlePrimary);
        const pandora::PfoList &parentList(pPfo->GetParentPfoList());

        if (!parentList.empty())
        {
            if (parentList.size() != 1)
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- this particle has multiple parent particles ";

            ThreeDParticleMap::const_iterator parentIdIter = particleMap.find(*parentList.begin());
            if (particleMap.end() == parentIdIter)
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- found an unassociated particle ";

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
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- found an unassociated particle ";

            const size_t daughterIdCode(daughterIdIter->second);
            daughterIdCodes.push_back(daughterIdCode);
        }

        // Build Particle
        recob::PFParticle newParticle(pPfo->GetParticleId(), pfoIdCode, parentIdCode, daughterIdCodes);
        outputParticles->push_back(newParticle);

        // Associate Vertex
        if (!pPfo->GetVertexList().empty())
        {
            if (pPfo->GetVertexList().size() != 1)
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- this particle has multiple interaction vertices ";

            const pandora::Vertex *const pVertex = *(pPfo->GetVertexList().begin());

            ThreeDVertexMap::const_iterator iter = vertexMap.find(pVertex);
            if (vertexMap.end() == iter)
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- found an unassociated vertex ";

            const unsigned int vtxElement(iter->second);
            util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputVertices.get()), *(outputParticlesToVertices.get()),
                vtxElement, vtxElement + 1);
        }

        // Build 2D Clusters
        pandora::ClusterVector pandoraClusterVector(pPfo->GetClusterList().begin(), pPfo->GetClusterList().end());
        std::sort(pandoraClusterVector.begin(), pandoraClusterVector.end(), lar_content::LArClusterHelper::SortByNHits);

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

                const geo::WireID wireID(hit->WireID());
                const unsigned int volID(100000 * wireID.Cryostat + wireID.TPC);
                hitArray[volID].push_back(hit);

                if (pCaloHit2D->IsIsolated())
                    isolatedHits.insert(hit);
            }

            if (hitArray.empty())
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- found a cluster with no hits ";

            for (const HitArray::value_type &hitArrayEntry : hitArray)
            {
                const HitVector &clusterHits(hitArrayEntry.second);
                outputClusters->emplace_back(LArPandoraOutput::BuildCluster(clusterCounter++, clusterHits, isolatedHits, ClusterParamAlgo));

                util::CreateAssn(*(settings.m_pProducer), evt, *(outputClusters.get()), clusterHits, *(outputClustersToHits.get()));
                util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputClusters.get()), *(outputParticlesToClusters.get()),
                    outputClusters->size() - 1, outputClusters->size());

                LOG_DEBUG("LArPandora") << "Stored cluster ID=" << outputClusters->back().ID() << " (#" << (outputClusters->size() - 1)
                    << ") with " << clusterHits.size() << " hits";
            }
        }

        // Build 3D SpacePoints and calculate T0 from shift in 2D hits
        pandora::CaloHitList pandoraHitList3D;
        lar_content::LArPfoHelper::GetCaloHits(pPfo, pandora::TPC_3D, pandoraHitList3D);

        pandora::CaloHitVector pandoraHitVector3D(pandoraHitList3D.begin(), pandoraHitList3D.end());
        std::sort(pandoraHitVector3D.begin(), pandoraHitVector3D.end(), lar_content::LArClusterHelper::SortHitsByPosition);

        double sumT(0.), sumN(0.);

        for (const pandora::CaloHit *const pCaloHit3D : pandoraHitVector3D)
        {
            if (pandora::TPC_3D != pCaloHit3D->GetHitType())
                throw cet::exception("LArPandora") << " LArPandoraOutput::ProduceArtOutput --- found a 2D hit in a 3D cluster";

            const pandora::CaloHit *const pCaloHit2D = static_cast<const pandora::CaloHit*>(pCaloHit3D->GetParentAddress());

            const art::Ptr<recob::Hit> hit = LArPandoraOutput::GetHit(idToHitMap, pCaloHit2D);

            HitVector spacePointHits;
            spacePointHits.push_back(hit);

            // ATTN: We assume that the 2D Pandora hits have been shifted
            sumT += LArPandoraOutput::CalculateT0(hit, pCaloHit2D);
            sumN += 1.;

            outputSpacePoints->emplace_back(LArPandoraOutput::BuildSpacePoint(spacePointCounter++, pCaloHit3D));

            util::CreateAssn(*(settings.m_pProducer), evt, *(outputSpacePoints.get()), spacePointHits, *(outputSpacePointsToHits.get()));
            util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputSpacePoints.get()), *(outputParticlesToSpacePoints.get()),
                outputSpacePoints->size() - 1, outputSpacePoints->size());
        }

        // Output T0 objects [arguments are:  time (nanoseconds);  trigger type (3 for TPC stitching!);  pfparticle SelfID code;  T0 ID code]
        // ATTN: T0 values are currently calculated in nanoseconds relative to the trigger offset. Only non-zero values are outputted.
        const double T0((sumN > 0. && std::fabs(sumT) > sumN) ? (sumT / sumN) : 0.);

        if (settings.m_buildStitchedParticles && std::fabs(T0) > 0.)
        {
            outputT0s->emplace_back(anab::T0(T0, 3, outputParticles->back().Self(), t0Counter++));
            util::CreateAssn(*(settings.m_pProducer), evt, *(outputParticles.get()), *(outputT0s.get()), *(outputParticlesToT0s.get()), outputT0s->size() - 1, outputT0s->size());
        }
    }

    mf::LogDebug("LArPandora") << "   Number of new particles: " << outputParticles->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new clusters: " << outputClusters->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new space points: " << outputSpacePoints->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new vertices: " << outputVertices->size() << std::endl;

    if (settings.m_buildStitchedParticles)
        mf::LogDebug("LArPandora") << "   Number of new T0s: " << outputT0s->size() << std::endl;

    evt.put(std::move(outputParticles));
    evt.put(std::move(outputSpacePoints));
    evt.put(std::move(outputClusters));
    evt.put(std::move(outputVertices));

    evt.put(std::move(outputParticlesToSpacePoints));
    evt.put(std::move(outputParticlesToClusters));
    evt.put(std::move(outputParticlesToVertices));
    evt.put(std::move(outputSpacePointsToHits));
    evt.put(std::move(outputClustersToHits));

    if (settings.m_buildStitchedParticles)
    {
        evt.put(std::move(outputT0s));
        evt.put(std::move(outputParticlesToT0s));
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

recob::SpacePoint LArPandoraOutput::BuildSpacePoint(const int id, const pandora::CaloHit *const pCaloHit)
{
    if (pandora::TPC_3D != pCaloHit->GetHitType())
        throw cet::exception("LArPandora") << " LArPandoraOutput::BuildSpacePoint --- trying to build a space point from a 2D hit";

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
        throw cet::exception("LArPandora") << " LArPandoraOutput::GetHit --- found a Pandora hit without a parent ART hit ";

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

LArPandoraOutput::Settings::Settings() :
    m_pPrimaryPandora(nullptr),
    m_pProducer(nullptr),
    m_buildStitchedParticles(false)
{
}

} // namespace lar_pandora
