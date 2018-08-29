/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraExternalEventBuilding.cc
 *
 *  @brief  module for lar pandora external event building
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larpandora/LArPandoraEventBuilding/Slice.h"
#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraEvent.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "TTree.h"

namespace lar_pandora
{

class LArPandoraExternalEventBuilding : public art::EDProducer
{
public:
    explicit LArPandoraExternalEventBuilding(fhicl::ParameterSet const & pset);
    
    LArPandoraExternalEventBuilding(LArPandoraExternalEventBuilding const &) = delete;
    LArPandoraExternalEventBuilding(LArPandoraExternalEventBuilding &&) = delete;
    LArPandoraExternalEventBuilding & operator = (LArPandoraExternalEventBuilding const &) = delete;
    LArPandoraExternalEventBuilding & operator = (LArPandoraExternalEventBuilding &&) = delete;

    void produce(art::Event &evt) override;
    void endSubRun(art::SubRun &subrun);

private:
    typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<larpandoraobj::PFParticleMetadata> > PFParticleToMetadata;

    /**
     *  @brief  Collect PFParticles from the ART event and their mapping to metadata objects
     *
     *  @param  evt the ART event
     *  @param  particlesToMetadata the output mapping from PFParticles to their metadata
     *  @param  particles the output vector of particles
     */
    void CollectPFParticles(const art::Event &evt, PFParticleToMetadata &particlesToMetadata, PFParticleVector &particles) const;

    /**
     *  @brief  Build mapping from ID to PFParticle for fast navigation through the hierarchy
     *
     *  @param  particlesToMetadata the input mapping from PFParticles to their metadata
     *  @param  particleMap the output mapping from ID to PFParticle
     */
    void BuildPFParticleMap(const PFParticleToMetadata &particlesToMetadata, PFParticleMap &particleMap) const;

    /**
     *  @brief  Collect PFParticles that have been identified as clear cosmic ray muons by pandora
     *
     *  @param  allParticles input vector of all particles
     *  @param  particlesToMetadata the input mapping from PFParticles to their metadata
     *  @param  particleMap the input mapping from ID to PFParticle
     *  @param  clearCosmics the output vector of clear cosmic rays
     */
    void CollectClearCosmicRays(const PFParticleVector &allParticles, const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, PFParticleVector &clearCosmics) const;

    /**
     *  @brief  Collect slices 
     *
     *  @param  allParticles input vector of all particles
     *  @param  particlesToMetadata the input mapping from PFParticles to their metadata
     *  @param  particleMap the input mapping from ID to PFParticle
     *  @param  slices the output vector of slices
     */
    void CollectSlices(const PFParticleVector &allParticles, const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, SliceVector &slices) const;

    /**
     *  @brief  Get the consolidated collection of particles based on the slice ids
     *
     *  @param  allParticles input vector of all particles
     *  @param  clearCosmics the input vector of clear cosmic ray muons
     *  @param  slices the input vector of slices
     *  @param  consolidatedParticles the output vector of particles to include in the consolidated output 
     */
    void CollectConsolidatedParticles(const PFParticleVector &allParticles, const PFParticleVector &clearCosmics, const SliceVector &slices, PFParticleVector &consolidatedParticles) const;

    /**
     *  @brief  Query a metadata object for a given key and return the corresponding value
     *
     *  @param  metadata the metadata object to query
     *  @param  key the key to search for
     *  
     *  @return the value in the metadata corresponding to the input key
     */
    float GetMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &key) const;

    /**
     *  @brief  Query a metadata object to see if it is a target particle
     *
     *  @param  metadata the metadata object to query
     *  
     *  @return boolean - if the particle is a target
     */
    bool IsTarget(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata) const;

    /**
     *  @brief  Query a metadata object to see if it is a target particle
     *
     *  @param  metadata the metadata object to query
     *  
     *  @return boolean - if the particle is a target
     */
    bool IsTarget(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata) const;

    std::string                         m_inputProducerLabel;  ///< Label for the Pandora instance that produced the collections we want to consolidated
    std::string                         m_trackProducerLabel;  ///< Label for the track producer using the Pandora instance that produced the collections we want to consolidate
    std::string                         m_showerProducerLabel; ///< Label for the shower producer using the Pandora instance that produced the collections we want to consolidate
    std::string                         m_hitProducerLabel;    ///< Label for the hit producer that was used as input to the Pandora instance specified
    bool                                m_shouldProduceT0s;    ///< If we should produce T0s (relevant when stitching over multiple drift volumes)
    art::InputTag                       m_pandoraTag;          ///< The input tag for the pandora producer
    std::unique_ptr<SliceIdBaseTool>    m_sliceIdTool;         ///< The slice id tool
    bool                                m_useTestBeamMode;     ///< If we should expect a test-beam (instead of a neutrino) slice
    std::string                         m_targetKey;           ///< The metadata key for a PFParticle to determine if it is the target
    std::string                         m_scoreKey;            ///< The metadata key for the score of the target slice from Pandora
    bool                                m_isData;              ///< If this is a data event
    std::string                         m_generatorLabel;      ///< The label of the generator for MC event POT counting
    int                                 m_run;                 ///< The run number
    int                                 m_subRun;              ///< The subRun number
    float                               m_pot;                 ///< The total amount of POT for the current sub run
    TTree                              *m_pSubRunTree;         ///< The tree holding subrun information for POT counting of MC samples
};

DEFINE_ART_MODULE(LArPandoraExternalEventBuilding)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "Pandora/PdgTable.h"
#include "larcoreobj/SummaryData/POTSummary.h"

namespace lar_pandora
{

LArPandoraExternalEventBuilding::LArPandoraExternalEventBuilding(fhicl::ParameterSet const &pset) :
    m_inputProducerLabel(pset.get<std::string>("InputProducerLabel")),
    m_trackProducerLabel(pset.get<std::string>("TrackProducerLabel")),
    m_showerProducerLabel(pset.get<std::string>("ShowerProducerLabel")),
    m_hitProducerLabel(pset.get<std::string>("HitProducerLabel")),
    m_shouldProduceT0s(pset.get<bool>("ShouldProduceT0s")),
    m_pandoraTag(art::InputTag(m_inputProducerLabel)),
    m_sliceIdTool(art::make_tool<SliceIdBaseTool>(pset.get<fhicl::ParameterSet>("SliceIdTool"))),
    m_useTestBeamMode(pset.get<bool>("ShouldUseTestBeamMode", false)),
    m_targetKey(m_useTestBeamMode ? "IsTestBeam" : "IsNeutrino"),
    m_scoreKey(m_useTestBeamMode ? "TestBeamScore" : "NuScore"),
    m_isData(pset.get<bool>("IsData")),
    m_generatorLabel(m_isData ? "" : pset.get<std::string>("GeneratorLabel")),
    m_run(std::numeric_limits<unsigned int>::max()),
    m_subRun(std::numeric_limits<unsigned int>::max()),
    m_pot(-std::numeric_limits<float>::max()),
    m_pSubRunTree(nullptr)
{
    produces< std::vector<recob::PFParticle> >();
    produces< std::vector<recob::SpacePoint> >();
    produces< std::vector<recob::Cluster> >();
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::Track> >(); 
    produces< std::vector<recob::Shower> >();
    produces< std::vector<recob::PCAxis> >();
    produces< std::vector<larpandoraobj::PFParticleMetadata> >();

    produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
    produces< art::Assns<recob::PFParticle, recob::Cluster> >();
    produces< art::Assns<recob::PFParticle, recob::Vertex> >();
    produces< art::Assns<recob::PFParticle, recob::Track> >();
    produces< art::Assns<recob::PFParticle, recob::Shower> >();
    produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
    produces< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();
    produces< art::Assns<recob::Shower, recob::PCAxis> >();
    produces< art::Assns<recob::SpacePoint, recob::Hit> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();

    if (m_shouldProduceT0s)
    {
        produces< std::vector<anab::T0> >();
        produces< art::Assns<recob::PFParticle, anab::T0> >();
    }

    if (!m_isData)
    {
        art::ServiceHandle<art::TFileService> fileService;

        m_pSubRunTree = fileService->make<TTree>("subruns","");
        m_pSubRunTree->Branch("run"   , &m_run   , "run/I");
        m_pSubRunTree->Branch("subRun", &m_subRun, "subRun/I");
        m_pSubRunTree->Branch("pot"   , &m_pot   , "pot/F");
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::produce(art::Event &evt)
{
    PFParticleVector particles;
    PFParticleToMetadata particlesToMetadata;
    this->CollectPFParticles(evt, particlesToMetadata, particles);

    PFParticleMap particleMap;
    this->BuildPFParticleMap(particlesToMetadata, particleMap);

    PFParticleVector clearCosmics;
    this->CollectClearCosmicRays(particles, particlesToMetadata, particleMap, clearCosmics);

    SliceVector slices;
    this->CollectSlices(particles, particlesToMetadata, particleMap, slices);
    
    m_sliceIdTool->ClassifySlices(slices, evt);

    PFParticleVector consolidatedParticles;
    this->CollectConsolidatedParticles(particles, clearCosmics, slices, consolidatedParticles);

    const LArPandoraEvent::Labels labels(m_inputProducerLabel, m_trackProducerLabel, m_showerProducerLabel, m_hitProducerLabel); 
    const LArPandoraEvent consolidatedEvent(LArPandoraEvent(this, &evt, labels, m_shouldProduceT0s), consolidatedParticles);

    consolidatedEvent.WriteToEvent();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::CollectPFParticles(const art::Event &evt, PFParticleToMetadata &particlesToMetadata, PFParticleVector &particles) const
{
    art::Handle<std::vector<recob::PFParticle> > pfParticleHandle;
    evt.getByLabel(m_pandoraTag, pfParticleHandle);

    art::FindManyP<larpandoraobj::PFParticleMetadata> pfParticleMetadataAssoc(pfParticleHandle, evt, m_pandoraTag);
  
    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> part(pfParticleHandle, i);
        const auto &metadata(pfParticleMetadataAssoc.at(part.key()));

        particles.push_back(part);

        if (metadata.size() != 1) 
            throw cet::exception("LArPandora") << " LArPandoraExternalEventBuilding::CollectPFParticles -- Found a PFParticle without exactly 1 metadata associated." << std::endl;

        if (!particlesToMetadata.insert(PFParticleToMetadata::value_type(part, metadata.front())).second)
            throw cet::exception("LArPandoraExternalEventBuilding") << "Repeated PFParticles" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::BuildPFParticleMap(const PFParticleToMetadata &particlesToMetadata, PFParticleMap &particleMap) const
{
    for (const auto &entry : particlesToMetadata)
    {
        if (!particleMap.insert(PFParticleMap::value_type(entry.first->Self(), entry.first)).second)
            throw cet::exception("LArPandoraExternalEventBuilding") << "Repeated PFParticles" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::CollectClearCosmicRays(const PFParticleVector &allParticles, const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, PFParticleVector &clearCosmics) const
{
    for (const auto &part : allParticles)
    {
        // Get the parent of the particle
        const auto parentIt(particlesToMetadata.find(LArPandoraHelper::GetParentPFParticle(particleMap, part)));
        if (parentIt == particlesToMetadata.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Found PFParticle without metadata" << std::endl;

        // ATTN particles without the "IsClearCosmic" parameter are not clear cosmics
        try
        {
            if (static_cast<bool>(std::round(this->GetMetadataValue(parentIt->second, "IsClearCosmic"))))
                clearCosmics.push_back(part);
        }
        catch (const cet::exception &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::CollectSlices(const PFParticleVector &allParticles, const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, SliceVector &slices) const
{
    std::map<unsigned int, float> targetScores;
    std::map<unsigned int, PFParticleVector> crHypotheses;
    std::map<unsigned int, PFParticleVector> targetHypotheses;
    std::vector<unsigned int> usedSliceIds;

    // Collect the slice information
    for (const auto &part : allParticles)
    {
        // Find the parent PFParticle
        const auto parentIt(particlesToMetadata.find(LArPandoraHelper::GetParentPFParticle(particleMap, part)));
        if (parentIt == particlesToMetadata.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Found PFParticle without metadata" << std::endl;
       
        // Skip PFParticles that are clear cosmics
        try
        {
            if (static_cast<bool>(std::round(this->GetMetadataValue(parentIt->second, "IsClearCosmic"))))
                continue;
        }
        catch (const cet::exception &)
        {
        }

        const unsigned int sliceId(static_cast<unsigned int>(std::round(this->GetMetadataValue(parentIt->second, "SliceIndex"))));
        const float targetScore(this->GetMetadataValue(parentIt->second, m_scoreKey));

        // Keep track of the slice IDs we have used, and their corresponding score
        if (std::find(usedSliceIds.begin(), usedSliceIds.end(), sliceId) == usedSliceIds.end())
        {
            usedSliceIds.push_back(sliceId);

            // ATTN all PFParticles in the same slice will have the same targetScore
            targetScores[sliceId] = targetScore;
        }

        if (this->IsTarget(parentIt->second))
        {
            targetHypotheses[sliceId].push_back(part);
        }
        else 
        {
            crHypotheses[sliceId].push_back(part);
        }
    }

    // Sort the slice IDs to ensure reproducibility
    std::sort(usedSliceIds.begin(), usedSliceIds.end());

    // ATTN: we need to ensure that for each slice there is a cosmic and neutrino hypothesis, even if the pass created no PFOs
    // in such a case we add an empty vector of pfparticles
    const PFParticleVector emptyPFParticleVector;

    // Produce the slices

    // ATTN slice indices are enumerated from 1
    for (unsigned int sliceId = 1; sliceId <= targetScores.size(); ++sliceId)
    {
        // Get the target score
        const auto targetScoresIter(targetScores.find(sliceId));
        if (targetScoresIter == targetScores.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Scrambled slice information - can't find target score with id = " << sliceId << std::endl;

        PFParticleVector targetPFParticleVector, crPFParticleVector;

        // Get the target hypothesis
        const auto targetHypothesisIter(targetHypotheses.find(sliceId));
        targetPFParticleVector = ((targetHypothesisIter == targetHypotheses.end()) ? emptyPFParticleVector : targetHypothesisIter->second);

        // Get the cosmic hypothesis
        const auto crHypothesisIter(crHypotheses.find(sliceId));
        crPFParticleVector = ((crHypothesisIter == crHypotheses.end()) ? emptyPFParticleVector : crHypothesisIter->second);
        slices.emplace_back(targetScoresIter->second, targetPFParticleVector, crPFParticleVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPandoraExternalEventBuilding::GetMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &key) const
{
    const auto &propertiesMap(metadata->GetPropertiesMap());
    const auto &it(propertiesMap.find(key));

    if (it == propertiesMap.end())
        throw cet::exception("LArPandoraExternalEventBuilding") << "No key \"" << key << "\" found in metadata properties map" << std::endl;

    return it->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::CollectConsolidatedParticles(const PFParticleVector &allParticles, const PFParticleVector &clearCosmics, const SliceVector &slices, PFParticleVector &consolidatedParticles) const
{
    PFParticleVector collectedParticles;
    collectedParticles.insert(collectedParticles.end(), clearCosmics.begin(), clearCosmics.end());

    for (const auto &slice : slices)
    {
        const PFParticleVector &particles(slice.IsTaggedAsTarget() ? slice.GetTargetHypothesis() : slice.GetCosmicRayHypothesis());
        collectedParticles.insert(collectedParticles.end(), particles.begin(), particles.end());
    }

    // ATTN the collected particles are the ones we want to output, but here we loop over all particles to ensure that the consolidated 
    // particles have the same ordering.
    for (const auto &part : allParticles)
    {
        if (std::find(collectedParticles.begin(), collectedParticles.end(), part) != collectedParticles.end())
            consolidatedParticles.push_back(part);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraExternalEventBuilding::IsTarget(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata) const
{
    try
    {
        return static_cast<bool>(std::round(this->GetMetadataValue(metadata, m_targetKey)));
    }
    catch (const cet::exception &)
    {
        return false;
    }
}

<<<<<<< HEAD
//------------------------------------------------------------------------------------------------------------------------------------------
    
void LArPandoraExternalEventBuilding::endSubRun(art::SubRun &subrun)
{
    if (m_isData)
        return;

    if (!m_pSubRunTree)
        throw cet::exception("LArPandora") << " LArPandoraExternalEventBuilding::endSubRun -- output tree not configured." << std::endl;

    art::Handle<sumdata::POTSummary> potSummaryHandle;

    m_run = subrun.run();
    m_subRun = subrun.subRun();
    m_pot = subrun.getByLabel(m_generatorLabel, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    
    m_pSubRunTree->Fill();  
}

=======
>>>>>>> Extended neutrino id tool to the general slice id tool - now works for protoDUNE
} // namespace lar_pandora

