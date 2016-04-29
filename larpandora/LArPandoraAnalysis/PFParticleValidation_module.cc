/**
 *  @file   larpandora/LArPandoraInterface/PFParticleMonitoring_module.cc
 *
 *  @brief  Analysis module for created particles
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

// ROOT includes
#include "TTree.h"

// Local LArPandora includes
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleValidation class
 */
class PFParticleValidation : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     PFParticleValidation(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~PFParticleValidation();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);
     void reconfigure(fhicl::ParameterSet const &pset);

private:
    typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
    typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;

    /**
     *  @brief Performing matching between true and reconstructed particles
     *
     *  @param recoParticlesToHits the mapping from reconstructed particles to hits
     *  @param trueParticlesToHits the mapping from true particles to hits
     *  @param hitsToTrueParticles the mapping from hits to true particles
     *  @param particleMatching the output matches between all reconstructed and true particles
     */
    void GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const MCParticlesToHits &trueParticlesToHits,
        const HitsToMCParticles &hitsToTrueParticles, ParticleMatchingMap &particleMatchingMap) const;

     std::string  m_hitfinderLabel;         ///<
     std::string  m_clusterLabel;           ///<
     std::string  m_particleLabel;          ///<
     std::string  m_geantModuleLabel;       ///<
};

DEFINE_ART_MODULE(PFParticleValidation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
/**
 *  @file   LArPandora/PFParticleValidation.cxx
 *
 *  @brief  Implementation of the lar pandora analysis producer.
 *
 *  $Log: $
 */

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Track.h"

// std includes
#include <iostream>

namespace lar_pandora
{

PFParticleValidation::PFParticleValidation(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleValidation::~PFParticleValidation()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleValidation::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_clusterLabel = pset.get<std::string>("ClusterModule","pandora");
    m_hitfinderLabel = pset.get<std::string>("HitFinderModule","gaushit");
    m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleValidation::beginJob()
{
    mf::LogDebug("LArPandora") << " *** PFParticleValidation::beginJob() *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleValidation::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleValidation::analyze(const art::Event &evt)
{
    // Collect Hits
    HitVector hitVector;
    LArPandoraHelper::CollectHits(evt, m_hitfinderLabel, hitVector);

    // Collect PFParticles and match Reco Particles to Hits
    PFParticlesToHits recoParticlesToHits;
    HitsToPFParticles recoHitsToParticles;
    LArPandoraHelper::BuildPFParticleHitMaps(evt, m_particleLabel, m_clusterLabel, recoParticlesToHits, recoHitsToParticles, LArPandoraHelper::kAddDaughters);

    // Collect MCParticles and match True Particles to Hits
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles hitsToTrueParticles;
    LArPandoraHelper::BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector, trueParticlesToHits, hitsToTrueParticles, LArPandoraHelper::kAddDaughters);

    // Match Reco Particles to True Particles
    ParticleMatchingMap particleMatchingMap;
    this->GetRecoToTrueMatches(recoParticlesToHits, trueParticlesToHits, hitsToTrueParticles, particleMatchingMap);

    for (const ParticleMatchingMap::value_type &trueToRecoEntry : particleMatchingMap)
    {
        std::cout << "MCPDG " << trueToRecoEntry.first->PdgCode()
                  << ", #MCHits " << trueParticlesToHits.at(trueToRecoEntry.first).size() <<  std::endl;

        for (const RecoParticleToNMatchedHits::value_type &recoToNMatchedHitsEntry : trueToRecoEntry.second)
        {
            std::cout << "--RecoPDG " << recoToNMatchedHitsEntry.first->PdgCode()
                      << ", #PfoHits " << recoParticlesToHits.at(recoToNMatchedHitsEntry.first).size()
                      << ", #MatchedHits " << recoToNMatchedHitsEntry.second << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleValidation::GetRecoToTrueMatches(const PFParticlesToHits &recoParticlesToHits, const MCParticlesToHits &trueParticlesToHits,
    const HitsToMCParticles &hitsToTrueParticles, ParticleMatchingMap &particleMatchingMap) const
{
    // Create a placeholder entry for all mc particles with >0 hits
    for (const MCParticlesToHits::value_type &trueParticleToHitsEntry : trueParticlesToHits)
    {
        if (!trueParticleToHitsEntry.second.empty())
            (void) particleMatchingMap.insert(ParticleMatchingMap::value_type(trueParticleToHitsEntry.first, RecoParticleToNMatchedHits()));
    }

    // Store true to reco matching details
    for (const PFParticlesToHits::value_type &recoParticleToHits : recoParticlesToHits)
    {
        const art::Ptr<recob::PFParticle> pRecoParticle(recoParticleToHits.first);
        const HitVector &hitVector(recoParticleToHits.second);

        for (const art::Ptr<recob::Hit> pHit : hitVector)
        {
            HitsToMCParticles::const_iterator trueParticleIter = hitsToTrueParticles.find(pHit);

            if (hitsToTrueParticles.end() == trueParticleIter)
                continue;

            const art::Ptr<simb::MCParticle> pTrueParticle = trueParticleIter->second;    
            particleMatchingMap[pTrueParticle][pRecoParticle]++;
        }
    }
}

} //namespace lar_pandora
