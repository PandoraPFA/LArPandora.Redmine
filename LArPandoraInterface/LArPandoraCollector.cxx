/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraCollector.cxx
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */

// Framework includes
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FindOneP.h"

// LArSoft includes
#include "SimpleTypesAndConstants/RawTypes.h" // raw::TDCtick_t
#include "Geometry/Geometry.h"
#include "Utilities/TimeService.h"

// Pandora includes
#include "Objects/ParticleFlowObject.h"

// Local LArPandora includes
#include "LArPandoraInterface/LArPandoraCollector.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

void LArPandoraCollector::CollectWires(const art::Event &evt, const std::string label, WireVector &wireVector)
{
    art::Handle< std::vector<recob::Wire> > theWires;
    evt.getByLabel(label, theWires);

    if (!theWires.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find wires... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theWires->size() << " Wires " << std::endl;
    }

    for (unsigned int i = 0; i < theWires->size(); ++i)
    {
        const art::Ptr<recob::Wire> wire(theWires, i);
        wireVector.push_back(wire);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectHits(const art::Event &evt, const std::string label, HitVector &hitVector)
{
    art::Handle< std::vector<recob::Hit> > theHits;
    evt.getByLabel(label, theHits);

    if (!theHits.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find hits... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theHits->size() << " Hits " << std::endl;
    }

    for (unsigned int i = 0; i < theHits->size(); ++i)
    {
        const art::Ptr<recob::Hit> hit(theHits, i);
        hitVector.push_back(hit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectPFParticles(const art::Event &evt, const std::string label, PFParticleVector &particleVector)
{
    art::Handle< std::vector<recob::PFParticle> > theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles " << std::endl;
    }

    for (unsigned int i = 0; i < theParticles->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> particle(theParticles, i);
        particleVector.push_back(particle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectSpacePoints(const art::Event &evt, const std::string label, SpacePointVector &spacePointVector, 
    SpacePointsToHits &spacePointsToHits)
{
    HitsToSpacePoints hitsToSpacePoints;
    return LArPandoraCollector::CollectSpacePoints(evt, label, spacePointVector, spacePointsToHits, hitsToSpacePoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectSpacePoints(const art::Event &evt, const std::string label, SpacePointVector &spacePointVector, 
    SpacePointsToHits &spacePointsToHits, HitsToSpacePoints &hitsToSpacePoints)
{
    art::Handle< std::vector<recob::SpacePoint> > theSpacePoints;
    evt.getByLabel(label, theSpacePoints);
 
    if (!theSpacePoints.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find spacepoints... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theSpacePoints->size() << " SpacePoints " << std::endl;
    }

    art::FindOneP<recob::Hit> theHitAssns(theSpacePoints, evt, label);
    for (unsigned int i = 0; i < theSpacePoints->size(); ++i)
    {
        const art::Ptr<recob::SpacePoint> spacepoint(theSpacePoints, i);
        spacePointVector.push_back(spacepoint);
        const art::Ptr<recob::Hit> hit = theHitAssns.at(i);
        spacePointsToHits[spacepoint] = hit;
        hitsToSpacePoints[hit] = spacepoint;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectClusters(const art::Event &evt, const std::string label, ClusterVector &clusterVector, 
    ClustersToHits &clustersToHits)
{
    art::Handle< std::vector<recob::Cluster> > theClusters;
    evt.getByLabel(label, theClusters);
 
    if (!theClusters.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find clusters... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theClusters->size() << " Clusters " << std::endl;
    }

    art::FindManyP<recob::Hit> theHitAssns(theClusters, evt, label);
    for (unsigned int i = 0; i < theClusters->size(); ++i)
    {
        const art::Ptr<recob::Cluster> cluster(theClusters, i);
        clusterVector.push_back(cluster);

        const std::vector< art::Ptr<recob::Hit> > hits = theHitAssns.at(i);
        for (unsigned int j=0; j<hits.size(); ++j)
        {
            const art::Ptr<recob::Hit> hit = hits.at(j);
            clustersToHits[cluster].push_back(hit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectPFParticles(const art::Event &evt, const std::string label, PFParticleVector &particleVector,
    PFParticlesToSpacePoints &particlesToSpacePoints)
{
    art::Handle< std::vector<recob::PFParticle> > theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles " << std::endl;
    }

    art::FindManyP<recob::SpacePoint> theSpacePointAssns(theParticles, evt, label);
    for (unsigned int i = 0; i < theParticles->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> particle(theParticles, i);
        particleVector.push_back(particle);

        const std::vector< art::Ptr<recob::SpacePoint> > spacepoints = theSpacePointAssns.at(i);
        for (unsigned int j=0; j<spacepoints.size(); ++j)
        {
            const art::Ptr<recob::SpacePoint> spacepoint = spacepoints.at(j);
            particlesToSpacePoints[particle].push_back(spacepoint);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectPFParticles(const art::Event &evt, const std::string label, PFParticleVector &particleVector,
    PFParticlesToClusters &particlesToClusters)
{
    art::Handle< std::vector<recob::PFParticle> > theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find particles... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " PFParticles " << std::endl;
    }

    art::FindManyP<recob::Cluster> theClusterAssns(theParticles, evt, label);
    for (unsigned int i = 0; i < theParticles->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> particle(theParticles, i);
        particleVector.push_back(particle);

        const std::vector< art::Ptr<recob::Cluster> > clusters = theClusterAssns.at(i);
        for (unsigned int j=0; j<clusters.size(); ++j)
        {
            const art::Ptr<recob::Cluster> cluster = clusters.at(j);
            particlesToClusters[particle].push_back(cluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectShowers(const art::Event &evt, const std::string label, ShowerVector &showerVector,
    PFParticlesToShowers &particlesToShowers)
{
    art::Handle< std::vector<recob::Shower> > theShowers;
    evt.getByLabel(label, theShowers);

    if (!theShowers.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find showers... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theShowers->size() << " Showers " << std::endl;
    }

    art::FindOneP<recob::PFParticle> theParticleAssns(theShowers, evt, label);
    for (unsigned int i = 0; i < theShowers->size(); ++i)
    {
        const art::Ptr<recob::Shower> shower(theShowers, i);
        showerVector.push_back(shower);
        const art::Ptr<recob::PFParticle> particle = theParticleAssns.at(i);
        particlesToShowers[particle].push_back(shower);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectTracks(const art::Event &evt, const std::string label, TrackVector &trackVector,
    PFParticlesToTracks &particlesToTracks)
{
    art::Handle< std::vector<recob::Track> > theTracks;
    evt.getByLabel(label, theTracks);

    if (!theTracks.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find tracks... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theTracks->size() << " Tracks " << std::endl;
    }

    art::FindOneP<recob::PFParticle> theParticleAssns(theTracks, evt, label);
    for (unsigned int i = 0; i < theTracks->size(); ++i)
    {
        const art::Ptr<recob::Track> track(theTracks, i);
        trackVector.push_back(track);
        const art::Ptr<recob::PFParticle> particle = theParticleAssns.at(i);
        particlesToTracks[particle].push_back(track);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectTracks(const art::Event &evt, const std::string label, TrackVector &trackVector, TracksToHits &tracksToHits)
{
    art::Handle< std::vector<recob::Track> > theTracks;
    evt.getByLabel(label, theTracks);

    if (!theTracks.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find tracks... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theTracks->size() << " Tracks " << std::endl;
    }

    art::FindManyP<recob::Hit> theHitAssns(theTracks, evt, label);
    for (unsigned int i = 0; i < theTracks->size(); ++i)
    {
        const art::Ptr<recob::Track> track(theTracks, i);
        trackVector.push_back(track);

        const std::vector< art::Ptr<recob::Hit> > hits = theHitAssns.at(i);
        for (unsigned int j=0; j<hits.size(); ++j)
        {
            const art::Ptr<recob::Hit> hit = hits.at(j);          
            tracksToHits[track].push_back(hit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectSeeds(const art::Event &evt, const std::string label, SeedVector &seedVector, 
    PFParticlesToSeeds &particlesToSeeds)
{
    art::Handle< std::vector<recob::Seed> > theSeeds;
    evt.getByLabel(label, theSeeds);

    if (!theSeeds.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find seeds... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theSeeds->size() << " Seeds " << std::endl;
    }

    art::FindManyP<recob::PFParticle> theSeedAssns(theSeeds, evt, label);
    for (unsigned int i = 0; i < theSeeds->size(); ++i)
    {
        const art::Ptr<recob::Seed> seed(theSeeds, i);
        seedVector.push_back(seed);

        const std::vector< art::Ptr<recob::PFParticle> > particles = theSeedAssns.at(i);
        for (unsigned int j=0; j<particles.size(); ++j)
        {
            const art::Ptr<recob::PFParticle> particle = particles.at(j);          
            particlesToSeeds[particle].push_back(seed);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectVertices(const art::Event &evt, const std::string label, VertexVector &vertexVector,
    PFParticlesToVertices &particlesToVertices)
{
    art::Handle< std::vector<recob::Vertex> > theVertices;
    evt.getByLabel(label, theVertices);

    if (!theVertices.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find vertices... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theVertices->size() << " Vertices " << std::endl;
    }

    art::FindManyP<recob::PFParticle> theVerticesAssns(theVertices, evt, label);
    for (unsigned int i = 0; i < theVertices->size(); ++i)
    {
        const art::Ptr<recob::Vertex> vertex(theVertices, i);
        vertexVector.push_back(vertex);

        const std::vector< art::Ptr<recob::PFParticle> > particles = theVerticesAssns.at(i);
        for (unsigned int j=0; j<particles.size(); ++j)
        {
            const art::Ptr<recob::PFParticle> particle = particles.at(j);          
            particlesToVertices[particle].push_back(vertex);
        }
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::BuildPFParticleHitMaps(const PFParticleVector &particleVector, const PFParticlesToSpacePoints &particlesToSpacePoints, 
    const SpacePointsToHits &spacePointsToHits, PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, 
    const DaughterMode daughterMode)
{ 
    // Build mapping from particle to particle ID for parent/daughter navigation
    PFParticleMap particleMap;

    for (PFParticleVector::const_iterator iter1 = particleVector.begin(), iterEnd1 = particleVector.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = *iter1;
        particleMap[particle->Self()] = particle;
    }

    // Loop over hits and build mapping between reconstructed final-state particles and reconstructed hits
    for (PFParticlesToSpacePoints::const_iterator iter1 = particlesToSpacePoints.begin(), iterEnd1 = particlesToSpacePoints.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> thisParticle = iter1->first;
        const art::Ptr<recob::PFParticle> particle((kAddDaughters == daughterMode) ? 
            LArPandoraCollector::GetParentPFParticle(particleMap, thisParticle) : thisParticle);

        if ((kIgnoreDaughters == daughterMode) && !LArPandoraCollector::IsFinalState(particleMap, particle))
            continue;

        const SpacePointVector &spacePointVector = iter1->second;

        for (SpacePointVector::const_iterator iter2 = spacePointVector.begin(), iterEnd2 = spacePointVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::SpacePoint> spacepoint = *iter2;

            SpacePointsToHits::const_iterator iter3 = spacePointsToHits.find(spacepoint);
            if (spacePointsToHits.end() == iter3)
                throw cet::exception("LArPandora") << " PandoraCollector::BuildPFParticleHitMaps --- Found a space point without an associated hit ";

            const art::Ptr<recob::Hit> hit = iter3->second;
        
            particlesToHits[particle].push_back(hit);
            hitsToParticles[hit] = particle;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::BuildPFParticleHitMaps(const PFParticleVector &particleVector, const PFParticlesToClusters &particlesToClusters, 
    const ClustersToHits &clustersToHits, PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, 
    const DaughterMode daughterMode)
{ 
    // Build mapping from particle to particle ID for parent/daughter navigation
    PFParticleMap particleMap;

    for (PFParticleVector::const_iterator iter1 = particleVector.begin(), iterEnd1 = particleVector.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> particle = *iter1;
        particleMap[particle->Self()] = particle;
    }

    // Loop over hits and build mapping between reconstructed final-state particles and reconstructed hits
    for (PFParticlesToClusters::const_iterator iter1 = particlesToClusters.begin(), iterEnd1 = particlesToClusters.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> thisParticle = iter1->first;
        const art::Ptr<recob::PFParticle> particle((kAddDaughters == daughterMode) ? 
            LArPandoraCollector::GetParentPFParticle(particleMap, thisParticle) : thisParticle);

        if ((kIgnoreDaughters == daughterMode) && !LArPandoraCollector::IsFinalState(particleMap, particle))
            continue;

        const ClusterVector &clusterVector = iter1->second;
        for (ClusterVector::const_iterator iter2 = clusterVector.begin(), iterEnd2 = clusterVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Cluster> cluster = *iter2;

            ClustersToHits::const_iterator iter3 = clustersToHits.find(cluster);
            if (clustersToHits.end() == iter3)
                throw cet::exception("LArPandora") << " PandoraCollector::BuildPFParticleHitMaps --- Found a space point without an associated hit ";

            const HitVector &hitVector = iter3->second;
            for (HitVector::const_iterator iter4 = hitVector.begin(), iterEnd4 = hitVector.end(); iter4 != iterEnd4; ++iter4)
            { 
                const art::Ptr<recob::Hit> hit = *iter4;
        
                particlesToHits[particle].push_back(hit);
                hitsToParticles[hit] = particle;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::BuildPFParticleHitMaps(const art::Event &evt, const std::string label_pfpart, const std::string label_middle,
    PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, const DaughterMode daughterMode, const bool useClusters)
{
    // Use intermediate clusters
    if (useClusters)
    {    
        PFParticleVector particleVector;
        PFParticlesToClusters particlesToClusters;

        ClusterVector clusterVector; 
        ClustersToHits clustersToHits;

        LArPandoraCollector::CollectPFParticles(evt, label_pfpart, particleVector, particlesToClusters);
        LArPandoraCollector::CollectClusters(evt, label_middle, clusterVector, clustersToHits);

        LArPandoraCollector::BuildPFParticleHitMaps(particleVector, particlesToClusters, clustersToHits, 
            particlesToHits, hitsToParticles, daughterMode);
    }

    // Use intermediate space points
    else
    {
        PFParticleVector particleVector;
        PFParticlesToSpacePoints particlesToSpacePoints;

        SpacePointVector spacePointVector; 
        SpacePointsToHits spacePointsToHits;

        LArPandoraCollector::CollectPFParticles(evt, label_pfpart, particleVector, particlesToSpacePoints);
        LArPandoraCollector::CollectSpacePoints(evt, label_middle, spacePointVector, spacePointsToHits);

        LArPandoraCollector::BuildPFParticleHitMaps(particleVector, particlesToSpacePoints, spacePointsToHits, 
            particlesToHits, hitsToParticles, daughterMode);   
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::SelectFinalStatePFParticles(const PFParticleVector &inputParticles, PFParticleVector &outputParticles)
{
    // Build mapping from particle to particle ID for parent/daughter navigation
    PFParticleMap particleMap;

    for (PFParticleVector::const_iterator iter = inputParticles.begin(), iterEnd = inputParticles.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;
        particleMap[particle->Self()] = particle;
    }
   
    // Select final-state particles
    for (PFParticleVector::const_iterator iter = inputParticles.begin(), iterEnd = inputParticles.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::PFParticle> particle = *iter;

        if (LArPandoraCollector::IsFinalState(particleMap, particle))
            outputParticles.push_back(particle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectSimChannels(const art::Event &evt, const std::string label, SimChannelVector &simChannelVector)
{
    if (evt.isRealData())
        throw cet::exception("LArPandora") << " PandoraCollector::CollectSimChannels --- Trying to access MC truth from real data ";

    art::Handle< std::vector<sim::SimChannel> > theSimChannels;
    evt.getByLabel(label, theSimChannels);
    
    if (!theSimChannels.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find sim channels... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theSimChannels->size() << " SimChannels " << std::endl;
    }

    for (unsigned int i = 0; i < theSimChannels->size(); ++i)
    {
        const art::Ptr<sim::SimChannel> channel(theSimChannels, i);
        simChannelVector.push_back(channel);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectMCParticles(const art::Event &evt, const std::string label, MCParticleVector &particleVector)
{
    if (evt.isRealData())
        throw cet::exception("LArPandora") << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

    art::Handle< std::vector<simb::MCParticle> > theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
    }

    for (unsigned int i = 0; i < theParticles->size(); ++i)
    {
        const art::Ptr<simb::MCParticle> particle(theParticles, i);
        particleVector.push_back(particle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::CollectMCParticles(const art::Event &evt, const std::string label, MCTruthToMCParticles &truthToParticles,
    MCParticlesToMCTruth &particlesToTruth)
{
    if (evt.isRealData())
        throw cet::exception("LArPandora") << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

    art::Handle< std::vector<simb::MCParticle> > theParticles;
    evt.getByLabel(label, theParticles);

    if (!theParticles.isValid())
    {
        mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
        return;
    }
    else
    {
        mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
    }

    art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);

    for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
    {
        const art::Ptr<simb::MCParticle> particle(theParticles, i);
        const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
        truthToParticles[truth].push_back(particle);
        particlesToTruth[particle] = truth;
    }
}
 
//------------------------------------------------------------------------------------------------------------------------------------------


void LArPandoraCollector::BuildMCParticleHitMaps(const HitVector &hitVector, const SimChannelVector &simChannelVector, 
    HitsToTrackIDEs &hitsToTrackIDEs)
{
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::TimeService> ts;

    SimChannelVector sortedSimChannelVector(geom->Nchannels());

    for (SimChannelVector::const_iterator iter = simChannelVector.begin(), iterEnd = simChannelVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<sim::SimChannel> simChannel = *iter;
        sortedSimChannelVector.at(simChannel->Channel()) = simChannel;
    }

    for (HitVector::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;
        
        const raw::TDCtick_t start_tdc(ts->TPCTick2TDC(hit->PeakTimeMinusRMS()));
        const raw::TDCtick_t end_tdc(ts->TPCTick2TDC(hit->PeakTimePlusRMS()));

        const TrackIDEVector trackCollection(sortedSimChannelVector.at(hit->Channel())->TrackIDEs(start_tdc, end_tdc));

        if (trackCollection.empty())
            continue; // Hit has no truth information [continue]

        for (unsigned int iTrack = 0, iTrackEnd = trackCollection.size(); iTrack < iTrackEnd; ++iTrack)
        {
            const sim::TrackIDE trackIDE = trackCollection.at(iTrack);
            hitsToTrackIDEs[hit].push_back(trackIDE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::BuildMCParticleHitMaps(const HitsToTrackIDEs &hitsToTrackIDEs, const MCTruthToMCParticles &truthToParticles,
    MCParticlesToHits &particlesToHits, HitsToMCParticles &hitsToParticles, const DaughterMode daughterMode)
{
    // Build mapping between particles and track IDs for parent/daughter navigation
    MCParticleMap particleMap;

    for (MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(), iterEnd1 = truthToParticles.end(); iter1 != iterEnd1; ++iter1)
    {
        const MCParticleVector &particleVector = iter1->second;
        for (MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<simb::MCParticle> particle = *iter2;
            particleMap[particle->TrackId()] = particle;
        }
    }

    // Loop over hits and build mapping between reconstructed hits and true particles
    for (HitsToTrackIDEs::const_iterator iter1 = hitsToTrackIDEs.begin(), iterEnd1 = hitsToTrackIDEs.end(); iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::Hit> hit = iter1->first;
        const TrackIDEVector &trackCollection = iter1->second;

        int bestTrackID(-1);
        float bestEnergyFrac(0.f);

        for (TrackIDEVector::const_iterator iter2 = trackCollection.begin(), iterEnd2 = trackCollection.end(); iter2 != iterEnd2; ++iter2)
        {
            const sim::TrackIDE &trackIDE = *iter2;
            const int trackID(std::abs(trackIDE.trackID)); // TODO: Find out why std::abs is needed
            const float energyFrac(trackIDE.energyFrac);

            if (energyFrac > bestEnergyFrac)
            {
                bestEnergyFrac = energyFrac;
                bestTrackID = trackID;
            }
        }

        if (bestTrackID >= 0)
        {
            MCParticleMap::const_iterator iter3 = particleMap.find(bestTrackID);
            if (particleMap.end() == iter3)
                throw cet::exception("LArPandora") << " PandoraCollector::BuildMCParticleHitMaps --- Found a track ID without an MC Particle ";

            const art::Ptr<simb::MCParticle> thisParticle = iter3->second;
            const art::Ptr<simb::MCParticle> particle((kAddDaughters == daughterMode) ? 
                LArPandoraCollector::GetParentMCParticle(particleMap, thisParticle) : thisParticle);

            if ((kIgnoreDaughters == daughterMode) && !(particleMap.end() == particleMap.find(particle->Mother())))
                continue;

            particlesToHits[particle].push_back(hit);
            hitsToParticles[hit] = particle;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraCollector::BuildMCParticleHitMaps(const art::Event &evt, const std::string label, const HitVector &hitVector, 
    MCParticlesToHits &particlesToHits, HitsToMCParticles &hitsToParticles, const DaughterMode daughterMode)
{
    SimChannelVector simChannelVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    HitsToTrackIDEs hitsToTrackIDEs;

    LArPandoraCollector::CollectSimChannels(evt, label, simChannelVector);
    LArPandoraCollector::CollectMCParticles(evt, label, truthToParticles, particlesToTruth);
    LArPandoraCollector::BuildMCParticleHitMaps(hitVector, simChannelVector, hitsToTrackIDEs);
    LArPandoraCollector::BuildMCParticleHitMaps(hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> LArPandoraCollector::GetParentPFParticle(const PFParticleMap &particleMap, const art::Ptr<recob::PFParticle> inputParticle)
{
    int primaryTrackID(inputParticle->Self());

    if (!inputParticle->IsPrimary())
    {
        int parentTrackID(inputParticle->Parent());

        while(1)
        {
            PFParticleMap::const_iterator pIter1 = particleMap.find(parentTrackID); 
            if (particleMap.end() == pIter1)
                throw cet::exception("LArPandora") << " PandoraCollector::GetParentPFParticle --- Found a PFParticle with a particle ID ";

            const art::Ptr<recob::PFParticle> parentParticle = pIter1->second;
            if (LArPandoraCollector::IsNeutrino(parentParticle))
                break;

            primaryTrackID = parentTrackID;

            if (parentParticle->IsPrimary())
                break;

            parentTrackID = parentParticle->Parent();
        }
    }

    PFParticleMap::const_iterator pIter2 = particleMap.find(primaryTrackID);
    if (particleMap.end() == pIter2)
        throw cet::exception("LArPandora") << " PandoraCollector::GetParentPFParticle --- Found a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> outputParticle = pIter2->second;
    return outputParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<simb::MCParticle> LArPandoraCollector::GetParentMCParticle(const MCParticleMap &particleMap, const art::Ptr<simb::MCParticle> inputParticle)
{
    int primaryTrackID(inputParticle->TrackId());
    int parentTrackID(inputParticle->Mother());
   
    while(1)
    {
        MCParticleMap::const_iterator pIter1 = particleMap.find(parentTrackID); 
        if (particleMap.end() == pIter1)
            break; // Can't find MC Particle for this track ID [break]

        const art::Ptr<simb::MCParticle> particle = pIter1->second;

        primaryTrackID = parentTrackID;
        parentTrackID = particle->Mother();
    }

    MCParticleMap::const_iterator pIter2 = particleMap.find(primaryTrackID);
    if (particleMap.end() == pIter2)
        throw cet::exception("LArPandora") << " PandoraCollector::GetParentMCParticle --- Found a track ID without a MC particle ";

    const art::Ptr<simb::MCParticle> outputParticle = pIter2->second;
    return outputParticle;
}


//------------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Track> LArPandoraCollector::GetPrimaryTrack(const PFParticlesToTracks &particlesToTracks, const art::Ptr<recob::PFParticle> particle)
{
    PFParticlesToTracks::const_iterator tIter = particlesToTracks.find(particle);

    if (particlesToTracks.end() == tIter || tIter->second.empty())
        throw cet::exception("LArPandora") << " PandoraCollector::GetPrimaryTrack --- Failed to find associated track ";

    if (tIter->second.size() != 1)
        throw cet::exception("LArPandora") << " PandoraCollector::GetPrimaryTrack --- Found more than one associated track ";

    const art::Ptr<recob::Track> primaryTrack = *(tIter->second.begin());
    return primaryTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArPandoraCollector::GetParentNeutrino(const PFParticleMap &particleMap, const art::Ptr<recob::PFParticle> daughterParticle)
{
    art::Ptr<recob::PFParticle> parentParticle = LArPandoraCollector::GetParentPFParticle(particleMap, daughterParticle);

    if (LArPandoraCollector::IsNeutrino(parentParticle))
        return parentParticle->PdgCode();

    if (parentParticle->IsPrimary())
        return 0;

    const int parentID(parentParticle->Parent());

    PFParticleMap::const_iterator pIter = particleMap.find(parentID);
    if (particleMap.end() == pIter)
        throw cet::exception("LArPandora") << " PandoraCollector::GetParentNeutrino --- Found a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> neutrinoParticle = pIter->second;
    return neutrinoParticle->PdgCode();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraCollector::IsFinalState(const PFParticleMap &particleMap, const art::Ptr<recob::PFParticle> daughterParticle)
{
    if (LArPandoraCollector::IsNeutrino(daughterParticle))
        return false;

    if (daughterParticle->IsPrimary())
        return true;

    const int parentID(daughterParticle->Parent());

    PFParticleMap::const_iterator pIter = particleMap.find(parentID);
    if (particleMap.end() == pIter)
        throw cet::exception("LArPandora") << " PandoraCollector::IsFinalState --- Found a PFParticle without a particle ID ";

    const art::Ptr<recob::PFParticle> parentParticle = pIter->second;

    if (LArPandoraCollector::IsNeutrino(parentParticle))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraCollector::IsNeutrino(const art::Ptr<recob::PFParticle> particle)
{
    const int pdg(particle->PdgCode());

    // electron, muon, tau (use Pandora PDG tables)
    return ((pandora::NU_E == std::abs(pdg)) || (pandora::NU_MU == std::abs(pdg)) || (pandora::NU_TAU == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraCollector::IsTrack(const art::Ptr<recob::PFParticle> particle)
{
    const int pdg(particle->PdgCode());

    // muon, pion, proton, kaon (use Pandora PDG tables)
    return ((pandora::MU_MINUS == std::abs(pdg)) || (pandora::PI_PLUS == std::abs(pdg)) || (pandora::PROTON == std::abs(pdg)) || 
        (pandora::K_PLUS == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraCollector::IsShower(const art::Ptr<recob::PFParticle> particle)
{
    const int pdg(particle->PdgCode());

    // electron, photon (use Pandora PDG tables)
    return ((pandora::E_MINUS == std::abs(pdg)) || (pandora::PHOTON == std::abs(pdg)));
}

} // namespace lar_pandora
