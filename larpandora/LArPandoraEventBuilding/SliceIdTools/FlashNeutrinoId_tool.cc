/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "ubana/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
#include "larpandora/LArPandoraEventBuilding/Slice.h"

namespace lar_pandora
{

/**
 *  @brief  Neutrino ID tool that selects the most likely neutrino slice using PMT information
 */
class FlashNeutrinoId : SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  pset FHiCL parameter set
     */
    FlashNeutrinoId(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Classify slices as neutrino or cosmic
     *
     *  @param  slices the input vector of slices to classify
     *  @param  evt the art event
     */
    void ClassifySlices(SliceVector &slices, const art::Event &evt) override;

private:
    typedef std::vector<recob::OpFlash> FlashVector;

    /**
     *  @breif  Try to find the brightes flash with sufficent photoelectons that is in time with the beam
     *
     *  @param  evt the art event
     *  @param  beamFlash the output beam flash
     *
     *  @return if a suitable beam flash could be located
     */
    bool GetBeamFlash(const art::Event &evt, recob::OpFlash &beamFlash) const;

    /**
     *  @breif  Given an input list of flashes, determine which occured within the beam window
     *
     *  @param  flashes the input vector of all flashes
     *  @param  flashesInWindow the output vector of flashes that are in time with the beam
     */
    void GetFlashesInBeamWindow(const FlashVector &flashes, FlashVector &flashesInWindow) const;

    /**
     *  @breif  Get the 3D spacepoints (with charge) associated with the PFParticles in the slice that are produced from hits in the W view
     *
     *  @param  pfParticleMap the input mapping from PFParticle ID to PFParticle
     *  @param  pfParticleToSpacePointMap the input mapping from PFParticles to SpacePoints
     *  @param  spacePointToHitMap the input mapping from SpacePoints to Hits
     *  @param  slice the input slice
     *
     *  @return the output charged cluster
     */
    flashana::QCluster_t GetChargeCluster(const PFParticleMap &pfParticleMap, const PFParticlesToSpacePoints &pfParticleToSpacePointMap,
        const SpacePointsToHits &spacePointToHitMap, const Slice &slice) const;

    /**
     *  @breif  Collect all downstream particles of those in the input vector
     *
     *  @param  pfParticleMap the mapping from PFParticle ID to PFParticle
     *  @param  parentPFParticles the input vector of PFParticles
     *  @param  downstreamPFParticle the output vector of PFParticles including those downstream of the input
     */
    void CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const PFParticleVector &parentPFParticles, 
        PFParticleVector &downstreamPFParticles) const;

    /**
     *  @breif  Collect all downstream particles of a given particle
     *
     *  @param  pfParticleMap the mapping from PFParticle ID to PFParticle
     *  @param  particle the input PFParticle
     *  @param  downstreamPFParticle the output vector of PFParticles including those downstream of the input
     */
    void CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const art::Ptr<recob::PFParticle> &particle,
        PFParticleVector &downstreamPFParticles) const;

    std::string  m_flashLabel;       ///< The label of the flash producer
    std::string  m_pandoraLabel;     ///< The label of the allOutcomes pandora producer
    float        m_beamWindowStart;  ///< The start time of the beam window
    float        m_beamWindowEnd;    ///< The end time of the beam window
    float        m_minBeamFlashPE;   ///< The minimum number of photoelectrons required to consider a flash as the beam flash
};

DEFINE_ART_CLASS_TOOL(FlashNeutrinoId)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_pandora
{
    
FlashNeutrinoId::FlashNeutrinoId(fhicl::ParameterSet const &pset) :
    m_flashLabel(pset.get<std::string>("FlashLabel")),
    m_pandoraLabel(pset.get<std::string>("PandoraAllOutcomesLabel")),
    m_beamWindowStart(pset.get<float>("BeamWindowStartTime")),
    m_beamWindowEnd(pset.get<float>("BeamWindowEndTime")),
    m_minBeamFlashPE(pset.get<float>("BeamFlashPEThreshold"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::ClassifySlices(SliceVector &slices, const art::Event &evt) 
{
    if (slices.empty()) return;

    // Find the flash, if any, in time with the beam with the largest number of photoelectrons that is sufficiently bright
    recob::OpFlash beamFlash;
    if (!this->GetBeamFlash(evt, beamFlash))
        return;
 
    // Collect the PFParticles and their associations to SpacePoints
    PFParticleVector pfParticles;
    PFParticlesToSpacePoints pfParticleToSpacePointMap;
    LArPandoraHelper::CollectPFParticles(evt, m_pandoraLabel, pfParticles, pfParticleToSpacePointMap);

    // Collect the SpacePoints and their associations to Hits
    SpacePointVector spacePoints;
    SpacePointsToHits spacePointToHitMap;
    LArPandoraHelper::CollectSpacePoints(evt, m_pandoraLabel, spacePoints, spacePointToHitMap);

    // Build a map from PFParticle ID to PFParticle for navigation through the hierarchy
    PFParticleMap pfParticleMap;
    LArPandoraHelper::BuildPFParticleMap(pfParticles, pfParticleMap);

    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        const auto chargeCluster(this->GetChargeCluster(pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, slices.at(sliceIndex)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::GetBeamFlash(const art::Event &evt, recob::OpFlash &beamFlash) const
{
    // Collect all flashes from the event
    art::InputTag flashTag(m_flashLabel); 
    const auto flashes(*evt.getValidHandle<FlashVector >(flashTag));

    // Find those inside the beam window
    FlashVector flashesInBeamWindow;
    this->GetFlashesInBeamWindow(flashes, flashesInBeamWindow);

    if (flashesInBeamWindow.empty())
        return false;

    // Get the flash with the highest number of photoelectrons
    beamFlash = (*std::max_element(flashesInBeamWindow.begin(), flashesInBeamWindow.end(), [](const recob::OpFlash &a, const recob::OpFlash &b) {
        return a.TotalPE() < b.TotalPE();
    }));

    return (beamFlash.TotalPE() >= m_minBeamFlashPE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::GetFlashesInBeamWindow(const FlashVector &flashes, FlashVector &flashesInWindow) const
{
    for (const auto &flash : flashes)
    {
        const auto time(flash.Time());
        if (time > m_beamWindowStart && time < m_beamWindowEnd)
            flashesInWindow.push_back(flash);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::QCluster_t FlashNeutrinoId::GetChargeCluster(const PFParticleMap &pfParticleMap,
    const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap, const Slice &slice) const
{
    // Collect all PFParticles in the slice, including those downstream of the primaries
    // ATTN here we only use the neutrino hypothesis, in theory this should work with either (or indeed both with some thought)
    PFParticleVector allParticlesInSlice;
    this->CollectDownstreamPFParticles(pfParticleMap, slice.GetTargetHypothesis(), allParticlesInSlice);

    flashana::QCluster_t chargeCluster;
    for (const auto &particle : allParticlesInSlice)
    {
        // Get the associated spacepoints
        const auto &partToSpacePointIter(pfParticleToSpacePointMap.find(particle));
        if (partToSpacePointIter == pfParticleToSpacePointMap.end())
            continue;

        for (const auto &spacePoint : partToSpacePointIter->second)
        {
            // Get the associated hit
            const auto &spacePointToHitIter(spacePointToHitMap.find(spacePoint));
            if (spacePointToHitIter == spacePointToHitMap.end())
                continue;

            // Only use hits from the collection plane
            const auto &hit(spacePointToHitIter->second);
            if (hit->View() != geo::kZ)
                continue;
            
            // Add the charged point to the vector
            const auto &position(spacePoint->XYZ());
            chargeCluster.emplace_back(position[0], position[1], position[2], hit->Integral());
        }
    }

    return chargeCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const PFParticleVector &parentPFParticles,
    PFParticleVector &downstreamPFParticles) const
{
    for (const auto &particle : parentPFParticles)
        this->CollectDownstreamPFParticles(pfParticleMap, particle, downstreamPFParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::CollectDownstreamPFParticles(const PFParticleMap &pfParticleMap, const art::Ptr<recob::PFParticle> &particle,
    PFParticleVector &downstreamPFParticles) const
{
    if (std::find(downstreamPFParticles.begin(), downstreamPFParticles.end(), particle) == downstreamPFParticles.end())
        downstreamPFParticles.push_back(particle);

    for (const auto &daughterId : particle->Daughters())
    {
        const auto iter(pfParticleMap.find(daughterId));
        if (iter == pfParticleMap.end())
            throw cet::exception("FlashNeutrinoId") << "Scrambled PFParticle IDs" << std::endl;

        this->CollectDownstreamPFParticles(pfParticleMap, iter->second, downstreamPFParticles);
    }
}

} // namespace lar_pandora
