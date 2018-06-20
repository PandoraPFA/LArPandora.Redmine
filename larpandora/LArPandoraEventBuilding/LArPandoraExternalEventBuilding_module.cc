/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraExternalEventBuilding.cc
 *
 *  @brief  module for lar pandora external event building
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larpandora/LArPandoraEventBuilding/Slice.h"
#include "larpandora/LArPandoraEventBuilding/NeutrinoIdBaseTool.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

namespace lar_pandora
{

class LArPandoraExternalEventBuilding : public art::EDAnalyzer
{
public:
    explicit LArPandoraExternalEventBuilding(fhicl::ParameterSet const & pset);
    
    LArPandoraExternalEventBuilding(LArPandoraExternalEventBuilding const &) = delete;
    LArPandoraExternalEventBuilding(LArPandoraExternalEventBuilding &&) = delete;
    LArPandoraExternalEventBuilding & operator = (LArPandoraExternalEventBuilding const &) = delete;
    LArPandoraExternalEventBuilding & operator = (LArPandoraExternalEventBuilding &&) = delete;

    void analyze(art::Event const & evt) override;

private:
    typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<larpandoraobj::PFParticleMetadata> > PFParticleToMetadata;

    /**
     *  @brief  Collect PFParticles from the ART event and their mapping to metadata objects
     *
     *  @param  evt the ART event
     *  @param  particlesToMetadata the output mapping from PFParticles to their metadata
     */
    void CollectPFParticles(const art::Event &evt, PFParticleToMetadata &particlesToMetadata) const;

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
     *  @param  particlesToMetadata the input mapping from PFParticles to their metadata
     *  @param  particleMap the input mapping from ID to PFParticle
     *  @param  clearCosmics the output vector of clear cosmic rays
     */
    void CollectClearCosmicRays(const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, PFParticleVector &clearCosmics) const;

    /**
     *  @brief  Collect slices 
     *
     *  @param  particlesToMetadata the input mapping from PFParticles to their metadata
     *  @param  particleMap the input mapping from ID to PFParticle
     *  @param  slices the output vector of slices
     */
    void CollectSlices(const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, SliceVector &slices) const;

    /**
     *  @brief  Query a metadata object for a given key and return the corresponding value
     *
     *  @param  metadata the metadata object to query
     *  @param  key the key to search for
     *  
     *  @return the value in the metadata corresponding to the input key
     */
    float GetMetadataValue(const art::Ptr<larpandoraobj::PFParticleMetadata> &metadata, const std::string &key) const;

    art::InputTag                       m_pandoraTag;        ///< The input tag for the pandora producer
    std::unique_ptr<NeutrinoIdBaseTool> m_neutrinoIdTool;    ///< The neutrino id tool
};

DEFINE_ART_MODULE(LArPandoraExternalEventBuilding)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "Pandora/PdgTable.h"

namespace lar_pandora
{

LArPandoraExternalEventBuilding::LArPandoraExternalEventBuilding(fhicl::ParameterSet const &pset) :
    EDAnalyzer(pset),
    m_pandoraTag(art::InputTag(pset.get<std::string>("PandoraLabel"))),
    m_neutrinoIdTool(art::make_tool<NeutrinoIdBaseTool>(pset.get<fhicl::ParameterSet>("NeutrinoIdTool")))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::analyze(art::Event const &evt)
{
    PFParticleToMetadata particlesToMetadata;
    this->CollectPFParticles(evt, particlesToMetadata);

    PFParticleMap particleMap;
    this->BuildPFParticleMap(particlesToMetadata, particleMap);

    PFParticleVector clearCosmics;
    this->CollectClearCosmicRays(particlesToMetadata, particleMap, clearCosmics);

    SliceVector slices;
    this->CollectSlices(particlesToMetadata, particleMap, slices);

    // Output the information collected
    std::cout << "Found " << clearCosmics.size() << " clear cosmic rays" << std::endl;
    std::cout << "      " << slices.size() << " slices:" << std::endl;

    for (const auto &slice : slices)
    {
        std::cout << "        ---------------" << std::endl;
        std::cout << "        nuScore         = " << slice.GetNeutrinoScore() << std::endl;
        std::cout << "        nParticles (nu) = " << slice.GetNeutrinoHypothesis().size() << std::endl;
        std::cout << "        nParticles (cr) = " << slice.GetCosmicRayHypothesis().size() << std::endl;
    }

    m_neutrinoIdTool->ClassifySlices(slices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::CollectPFParticles(const art::Event &evt, PFParticleToMetadata &particlesToMetadata) const
{
    art::Handle<art::Assns<recob::PFParticle,larpandoraobj::PFParticleMetadata,void> > particleMetadataAssn;
    evt.getByLabel(m_pandoraTag, particleMetadataAssn);
   
    for (const auto &entry : *particleMetadataAssn)
    {
        if (!particlesToMetadata.insert(PFParticleToMetadata::value_type(entry.first, entry.second)).second)
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

void LArPandoraExternalEventBuilding::CollectClearCosmicRays(const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, PFParticleVector &clearCosmics) const
{
    for (const auto &entry : particlesToMetadata)
    {
        try
        {
            // Check if the parent PFParticle is clear cosmic
            const auto parentIterator(particlesToMetadata.find(LArPandoraHelper::GetParentPFParticle(particleMap, entry.first)));

            if (static_cast<bool>(std::round(this->GetMetadataValue(parentIterator->second, "IsClearCosmic"))))
                clearCosmics.push_back(entry.first);
        }
        catch (const cet::exception &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraExternalEventBuilding::CollectSlices(const PFParticleToMetadata &particlesToMetadata, const PFParticleMap &particleMap, SliceVector &slices) const
{
    std::map<unsigned int, float> nuScores;
    std::map<unsigned int, PFParticleVector> crHypotheses;
    std::map<unsigned int, PFParticleVector> nuHypotheses;

    // Collect the slice information
    for (const auto &entry : particlesToMetadata)
    {
        // Find the parent PFParticle
        const auto parentIterator(particlesToMetadata.find(LArPandoraHelper::GetParentPFParticle(particleMap, entry.first)));
        if (parentIterator == particlesToMetadata.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Can't find the parent of input PFParticle" << std::endl;
        
        const int pdg(std::abs(parentIterator->first->PdgCode())); 
        const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);

        unsigned int sliceId(std::numeric_limits<unsigned int>::max());
        float nuScore(-std::numeric_limits<float>::max()); 

        try
        {
            sliceId = static_cast<unsigned int>(std::round(this->GetMetadataValue(parentIterator->second, "SliceIndex")));
            nuScore = this->GetMetadataValue(parentIterator->second, "NuScore");
        }
        catch (const cet::exception &exception)
        {
            // The above information should only unavailable for clear cosmic PFParticles
            try
            {
                if (static_cast<bool>(std::round(this->GetMetadataValue(parentIterator->second, "IsClearCosmic")))) continue;
            }
            catch (const cet::exception &)
            {
            }

            throw exception;
        }

        // ATTN all PFParticles in the same slice will have the same nuScore
        nuScores[sliceId] = nuScore;

        if (isNeutrino)
        {
            nuHypotheses[sliceId].push_back(entry.first);
        }
        else 
        {
            crHypotheses[sliceId].push_back(entry.first);
        }
    }

    // Produce the slices
    for (unsigned int sliceId = 1; sliceId <= nuScores.size(); ++sliceId)
    {
        // Get the neutrino score
        const auto nuScoresIter(nuScores.find(sliceId));
        if (nuScoresIter == nuScores.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Scrambled slice information - can't find nuScore with id = " << sliceId << std::endl;
        
        // Get the neutrino hypothesis
        const auto nuHypothesisIter(nuHypotheses.find(sliceId));
        if (nuHypothesisIter == nuHypotheses.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Scrambled slice information - can't find neutrino hypothesis with id = " << sliceId << std::endl;
        
        // Get the cosmic hypothesis
        const auto crHypothesisIter(crHypotheses.find(sliceId));
        if (crHypothesisIter == crHypotheses.end())
            throw cet::exception("LArPandoraExternalEventBuilding") << "Scrambled slice information - can't find cosmic hypothesis with id = " << sliceId << std::endl;

        slices.emplace_back(nuScoresIter->second, nuHypothesisIter->second, crHypothesisIter->second);
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

} // namespace lar_pandora

