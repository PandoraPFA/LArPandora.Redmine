/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraCollector.h
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */
#ifndef LAR_PANDORA_COLLECTOR_H
#define LAR_PANDORA_COLLECTOR_H

// Framework includes
#include "art/Framework/Principal/Event.h"

// LArSoft includes
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/PFParticle.h"

#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "Simulation/SimChannel.h"
 
#include "AnalysisBase/CosmicTag.h"

namespace lar_pandora 
{

// Try to make my code more readable!
typedef std::vector< art::Ptr<recob::Wire> >        WireVector;
typedef std::vector< art::Ptr<recob::Hit> >         HitVector;
typedef std::vector< art::Ptr<recob::SpacePoint> >  SpacePointVector;
typedef std::vector< art::Ptr<recob::Cluster> >     ClusterVector;
typedef std::vector< art::Ptr<recob::Seed> >        SeedVector;
typedef std::vector< art::Ptr<recob::Vertex> >      VertexVector;
typedef std::vector< art::Ptr<recob::Track> >       TrackVector;
typedef std::vector< art::Ptr<recob::Shower> >      ShowerVector;
typedef std::vector< art::Ptr<recob::PFParticle> >  PFParticleVector;
typedef std::vector< art::Ptr<simb::MCTruth> >      MCTruthVector;
typedef std::vector< art::Ptr<simb::MCParticle> >   MCParticleVector;
typedef std::vector< art::Ptr<sim::SimChannel> >    SimChannelVector;
typedef std::vector< sim::TrackIDE >                TrackIDEVector;
typedef std::vector< art::Ptr<anab::CosmicTag> >    CosmicTagVector;

typedef std::map< art::Ptr<recob::PFParticle>, TrackVector >                  PFParticlesToTracks;
typedef std::map< art::Ptr<recob::PFParticle>, ShowerVector >                 PFParticlesToShowers;
typedef std::map< art::Ptr<recob::PFParticle>, ClusterVector >                PFParticlesToClusters;
typedef std::map< art::Ptr<recob::PFParticle>, SeedVector >                   PFParticlesToSeeds;
typedef std::map< art::Ptr<recob::PFParticle>, VertexVector >                 PFParticlesToVertices;
typedef std::map< art::Ptr<recob::PFParticle>, SpacePointVector >             PFParticlesToSpacePoints;
typedef std::map< art::Ptr<recob::PFParticle>, HitVector >                    PFParticlesToHits;
typedef std::map< art::Ptr<recob::Track>,      HitVector >                    TracksToHits;
typedef std::map< art::Ptr<recob::Cluster>,    HitVector >                    ClustersToHits;
typedef std::map< art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit> >         SpacePointsToHits;
typedef std::map< art::Ptr<simb::MCTruth>,     MCParticleVector >             MCTruthToMCParticles;
typedef std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> >      MCParticlesToMCTruth;
typedef std::map< art::Ptr<simb::MCParticle>,  HitVector >                    MCParticlesToHits;
typedef std::map< art::Ptr<simb::MCParticle>,  art::Ptr<recob::PFParticle> >  MCParticlesToPFParticles;
typedef std::map< art::Ptr<recob::Hit>,        art::Ptr<recob::SpacePoint> >  HitsToSpacePoints;
typedef std::map< art::Ptr<recob::Hit>,        art::Ptr<recob::PFParticle> >  HitsToPFParticles;
typedef std::map< art::Ptr<recob::Hit>,        art::Ptr<simb::MCParticle> >   HitsToMCParticles;
typedef std::map< art::Ptr<recob::Hit>,        TrackIDEVector >               HitsToTrackIDEs;
typedef std::map< art::Ptr<recob::Track>,      CosmicTagVector >              TracksToCosmicTags;

typedef std::map< int, art::Ptr<recob::PFParticle> >  PFParticleMap;
typedef std::map< int, art::Ptr<recob::Cluster> >     ClusterMap;
typedef std::map< int, art::Ptr<recob::SpacePoint> >  SpacePointMap;
typedef std::map< int, art::Ptr<recob::Hit> >         HitMap;
typedef std::map< int, art::Ptr<simb::MCParticle> >   MCParticleMap;

class LArPandoraCollector 
{
public:
  
    enum DaughterMode 
    {
        kIgnoreDaughters = 0,    // Only use parent particles
        kUseDaughters = 1,       // Use both parent and daughter partcles
        kAddDaughters = 2        // Absorb daughter particles into parent particles
    };

    /**
     *  @brief Collect the reconstructed wires from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the Wire list in the event
     *  @param wireVector the ouput vector of Wire objects
     */
    static void CollectWires(const art::Event &evt, const std::string label, WireVector &wireVector);

    /**
     *  @brief Collect the reconstructed Hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the Hit list in the event
     *  @param hitVector the ouput vector of Hit objects
     */
    static void CollectHits(const art::Event &evt, const std::string label, HitVector &hitVector);

    /**
     *  @brief Collect the reconstructed PFParticles from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     */
    static void CollectPFParticles(const art::Event &evt, const std::string label, PFParticleVector &particleVector);

    /**
     *  @brief Collect the reconstructed SpacePoints and associated hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the SpacePoint list in the event
     *  @param spacePointVector the output vector of SpacePoint objects
     *  @param spacePointsToHits the output map from SpacePoint to Hit objects
     */
    static void CollectSpacePoints(const art::Event &evt, const std::string label, SpacePointVector &spacePointVector, 
        SpacePointsToHits &spacePointsToHits);   

    /**
     *  @brief Collect the reconstructed SpacePoints and associated hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the SpacePoint list in the event
     *  @param spacePointVector the output vector of SpacePoint objects
     *  @param spacePointsToHits the output map from SpacePoint to Hit objects
     *  @param hitsToSpacePoints the output map from Hit to SpacePoint objects
     */
    static void CollectSpacePoints(const art::Event &evt, const std::string label, SpacePointVector &spacePointVector, 
        SpacePointsToHits &spacePointsToHits, HitsToSpacePoints &hitsToSpacePoints);   

    /**
     *  @brief Collect the reconstructed Clusters and associated hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the SpacePoint list in the event
     *  @param clusterVector the output vector of Cluster objects
     *  @param clustersToHits the output map from Cluster to Hit objects
     */
    static void CollectClusters(const art::Event &evt, const std::string label, ClusterVector &clusterVector, 
        ClustersToHits &clustersToHits);   

    /**
     *  @brief Collect the reconstructed PFParticles and associated SpacePoints from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     *  @param particlesToSpacePoints the output map from PFParticle to SpacePoint objects
     */
    static void CollectPFParticles(const art::Event &evt, const std::string label, PFParticleVector &particleVector,
        PFParticlesToSpacePoints &particlesToSpacePoints);  

    /**
     *  @brief Collect the reconstructed PFParticles and associated Clusters from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param particleVector the output vector of PFParticle objects
     *  @param particlesToClusters the output map from PFParticle to Cluster objects
     */
    static void CollectPFParticles(const art::Event &evt, const std::string label, PFParticleVector &particleVector,
        PFParticlesToClusters &particlesToClusters);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Showers from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param showerVector the output vector of Shower objects
     *  @param particlesToShowers the output map from PFParticle to Shower objects
     */
    static void CollectShowers(const art::Event &evt, const std::string label, ShowerVector &showerVector,
        PFParticlesToShowers &particlesToShowers);   

    /**
     *  @brief Collect the reconstructed PFParticles and associated Tracks from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param trackVector the output vector of Track objects
     *  @param particlesToTracks the output map from PFParticle to Track objects
     */
    static void CollectTracks(const art::Event &evt, const std::string label, TrackVector &trackVector,
        PFParticlesToTracks &particlesToTracks);   

    /**
     *  @brief Collect the reconstructed Tracks and associated Hits from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param trackVector the output vector of Track objects
     *  @param tracksToHits the output map from Track to Hit objects
     */
    static void CollectTracks(const art::Event &evt, const std::string label, TrackVector &trackVector,
        TracksToHits &tracksToHits);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Seeds from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param seedVector the output vector of Seed objects
     *  @param particlesToSeeds the output map from PFParticle to Seed objects
     */
    static void CollectSeeds(const art::Event &evt, const std::string label, SeedVector &seedVector,
        PFParticlesToSeeds &particlesToSeeds);

    /**
     *  @brief Collect the reconstructed PFParticles and associated Vertices from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the PFParticle list in the event
     *  @param vertexVector the output vector of Vertex objects
     *  @param particlesToVertices the output map from PFParticle to Vertex objects
     */
    static void CollectVertices(const art::Event &evt, const std::string label, VertexVector &vertexVector,
        PFParticlesToVertices &particlesToVertices);

    /**
     *  @brief Build mapping between PFParticles and Hits using PFParticle/SpacePoint/Hit maps
     *
     *  @param particleVector the input vector of PFParticle objects
     *  @param particlesToSpacePoints the input map from PFParticle to SpacePoint objects
     *  @param spacePointsToHits the input map from SpacePoint to Hit objects
     *  @param particlesToHits the output map from PFParticle to Hit objects
     *  @param hitsToParticles the output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildPFParticleHitMaps(const PFParticleVector &particleVector, const PFParticlesToSpacePoints &particlesToSpacePoints, 
        const SpacePointsToHits &spacePointsToHits, PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, 
        const DaughterMode daughterMode = kUseDaughters);   

    /**
     *  @brief Build mapping between PFParticles and Hits using PFParticle/Cluster/Hit maps
     *
     *  @param particleVector the input vector of PFParticle objects
     *  @param particlesToClusters the input map from PFParticle to Cluster objects
     *  @param clustersToHits the input map from Cluster to Hit objects
     *  @param particlesToHits the output map from PFParticle to Hit objects
     *  @param hitsToParticles the output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildPFParticleHitMaps(const PFParticleVector &particleVector, const PFParticlesToClusters &particlesToClusters, 
        const ClustersToHits &clustersToHits, PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, 
        const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief Build mapping between PFParticles and Hits starting from ART event record
     *
     *  @param evt the ART event record
     *  @param label_pfpart the label for the PFParticle list in the event
     *  @param label_space the label for the Intermediate list in the event
     *  @param particlesToHits output map from PFParticle to Hit objects
     *  @param hitsToParticles output map from Hit to PFParticle objects
     *  @param daughterMode treatment of daughter particles in construction of maps
     *  @param useClusters choice of intermediate object (true for Clusters, false for SpacePoints)
     */
    static void BuildPFParticleHitMaps(const art::Event &evt, const std::string label_pfpart, const std::string label_mid,
        PFParticlesToHits &particlesToHits, HitsToPFParticles &hitsToParticles, const DaughterMode daughterMode = kUseDaughters,
        const bool useClusters = true);

    /**
     *  @brief Collect a vector of cosmic tags from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the cosmic tag information in the event
     *  @param cosmicTagVector output vector of CosmicTag objects
     *  @param tracksToCosmicTags output map from tracks to cosmic tags
     */
    static void CollectCosmicTags(const art::Event &evt, const std::string label, CosmicTagVector &cosmicTagVector, 
        TracksToCosmicTags &tracksToCosmicTags);

    /**
     *  @brief Collect a vector of SimChannel objects from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param simChannelVector output vector of SimChannel objects
     */
    static void CollectSimChannels(const art::Event &evt, const std::string label, SimChannelVector &simChannelVector);

    /**
     *  @brief Collect a vector of MCParticle objects from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param particleVector the output vector of MCParticle objects
     */
    static void CollectMCParticles(const art::Event &evt, const std::string label, MCParticleVector &particleVector);

    /**
     *  @brief Collect truth information from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param truthToParticles output map from MCTruth to MCParticle objects
     *  @param particlesToTruth output map from MCParticle to MCTruth objects
     */
    static void CollectMCParticles(const art::Event &evt, const std::string label, MCTruthToMCParticles &truthToParticles, 
        MCParticlesToMCTruth &particlesToTruth);

    /**
     *  @brief Collect the links from reconstructed hits to their true energy deposits 
     *
     *  @param hitVector the input vector of reconstructed hits
     *  @param simChannelVector the input vector of SimChannels
     *  @param hitsToTrackIDEs the out map from hits to true energy deposits
     */
    static void BuildMCParticleHitMaps(const HitVector &hitVector, const SimChannelVector &simChannelVector, HitsToTrackIDEs &hitsToTrackIDEs);

    /**
     *  @brief Build mapping between Hits and MCParticles, starting from Hit/TrackIDE/MCParticle information
     *
     *  @param hitsToTrackIDEs the input map from hits to true energy deposites
     *  @param truthToParticles the input map of truth information
     *  @param particlesToHits the mapping between true particles and reconstructed hits
     *  @param hitsToParticles the mapping between reconstructed hits and true particles
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildMCParticleHitMaps(const HitsToTrackIDEs &hitsToTrackIDEs, const MCTruthToMCParticles &truthToParticles,
        MCParticlesToHits &particlesToHits, HitsToMCParticles &hitsToParticles, const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief Build mapping between Hits and MCParticles, starting from ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param hitVector the input vector of reconstructed hits
     *  @param particlesToHits the output mapping between true particles and reconstructed hits
     *  @param hitsToParticles the output mapping between reconstructed hits and true particles
     *  @param daughterMode treatment of daughter particles in construction of maps
     */
    static void BuildMCParticleHitMaps(const art::Event &evt, const std::string label, const HitVector &hitVector, 
        MCParticlesToHits &particlesToHits, HitsToMCParticles &hitsToParticles, const DaughterMode daughterMode = kUseDaughters);

    /**
     *  @brief Select final-state reconstructed particles from a list of all reconstructed particles
     *
     *  @param inputParticles the input vector of all particles (it has to be all of them!)
     *  @param outputParticles the output vector of final-state particles
     */
    static void SelectFinalStatePFParticles(const PFParticleVector &inputParticles, PFParticleVector &outputParticles);    

    /**
     *  @brief Return the parent final-state particle by navigating up the chain of parent/daughter associations
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return the output parent final-state particle
     */
    static art::Ptr<recob::PFParticle> GetParentPFParticle(const PFParticleMap &particleMap, const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Return the parent final-state particle by navigating up the chain of parent/daughter associations
     *
     *  @param particleMap the mapping between true particle and true track ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return the output parent final-state particle
     */
    static art::Ptr<simb::MCParticle> GetParentMCParticle(const MCParticleMap &particleMap, const art::Ptr<simb::MCParticle> daughterParticle);

    /**
     *  @brief Return the primary track associated with a PFParticle
     *
     *  @param particlesToTracks the mapping between particles and tracks
     *  @param particle  the input particle
     */
    static art::Ptr<recob::Track> GetPrimaryTrack(const PFParticlesToTracks &particlesToTracks, const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief Return the parent neutrino PDG code (or zero for cosmics) for a given reconstructed particle
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return the PDG code of the parent neutrinos (or zero for cosmics)
     */
    static int GetParentNeutrino(const PFParticleMap &particleMap, const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Determine whether a particle has been reconstructed as a final-state particle
     *
     *  @param particleMap the mapping between reconstructed particle and particle ID
     *  @param daughterParticle the input daughter particle
     *
     *  @return true/false
     */
    static bool IsFinalState(const PFParticleMap &particleMap, const art::Ptr<recob::PFParticle> daughterParticle);

    /**
     *  @brief Determine whether a particle has been reconstructed as a neutrino
     *
     *  @param particle the input particle
     *
     *  @return true/false
     */
    static bool IsNeutrino(const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief Determine whether a particle has been reconstructed as track-like
     *
     *  @param particle the input particle
     *
     *  @return true/false
     */
    static bool IsTrack(const art::Ptr<recob::PFParticle> particle);

    /**
     *  @brief Determine whether a particle has been reconstructed as shower-like
     *
     *  @param particle the input particle
     *
     *  @return true/false
     */
    static bool IsShower(const art::Ptr<recob::PFParticle> particle);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_COLLECTOR_H
