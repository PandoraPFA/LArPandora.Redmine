/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "ubana/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubana/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
#include "larpandora/LArPandoraEventBuilding/Slice.h"

#include "Objects/CartesianVector.h"

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
    /**
     *  @brief  Data to describe an amount of charge deposited in a given 3D position
     */
    class Deposition
    {
    public:
        /**
         *  @brief  Default constructor
         *
         *  @param  x the x-component of the charge position
         *  @param  y the z-component of the charge position
         *  @param  z the z-component of the charge position
         *  @param  charge the charge deposited
         *  @param  nPhotons the estimated numer of photons produced
         */
        Deposition(const float x, const float y, const float z, const float charge, const float nPhotons);

        float m_x;         ///< The x-component of the charge position
        float m_y;         ///< The z-component of the charge position
        float m_z;         ///< The z-component of the charge position
        float m_charge;    ///< The charge deposited
        float m_nPhotons;  ///< The estimated numer of photons produced
    };

    typedef std::vector<Deposition> DepositionVector;
    typedef std::vector<recob::OpFlash> FlashVector;

    /**
     *  @brief  Get the ordered vector of optical detector IDs, use the remapping provided in FHiCL if required
     *
     *  @param  pset FHiCL parameter set
     */
    void GetOrderedOpDetVector(fhicl::ParameterSet const &pset);

    /**
     *  @breif  Try to find the brightest flash with sufficent photoelectons that is in time with the beam
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
    DepositionVector GetDepositionVector(const PFParticleMap &pfParticleMap, const PFParticlesToSpacePoints &pfParticleToSpacePointMap,
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

    /**
     *  @brief  Convert from deposited charge to number of photons for a given particle
     *
     *  @param  charge the input charge
     *  @param  particle the input particle
     *
     *  @return the number of photons
     */
    float GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &particle) const;

    /**
     *  @breif  Determine if a given slice is compatible with the beam flash by applying pre-selection cuts
     *
     *  @param  depositionVector the charge cluster from the space points in the slice
     *  @param  beamFlash the beam flash
     *
     *  @return if the depositionVector is compatible with the beamFlash
     */
    bool IsSliceCompatibleWithBeamFlash(const DepositionVector &depositionVector, const recob::OpFlash &beamFlash) const;

    /**
     *  @brief  Get the centroid of the input charge cluster, weighted by charge
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the charge weighted centroid
     */
    pandora::CartesianVector GetChargeWeightedCenter(const DepositionVector &depositionVector) const;

    /**
     *  @brief  Get the total charge from an input charge cluster
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the total charge
     */
    float GetTotalCharge(const DepositionVector &depositionVector) const;

    /**
     *  @breif  Apply flash matching between the input charge cluster and beam flash and return the score
     *
     *  @param  depositionVector the input charge cluster
     *  @param  beamFlash the input beam flash
     *
     *  @return the flash match score
     */
    float GetFlashMatchScore(const DepositionVector &depositionVector, const recob::OpFlash &beamFlash);

    /**
     *  @breif  Convert a recob::OpFlash into a flashana::Flash_t
     *
     *  @param  the input recob::OpFlash
     *
     *  @return the flashana::Flash_t
     */
    flashana::Flash_t ConvertFlashFormat(const recob::OpFlash &beamFlash) const;

    /**
     *  @brief  Convert a charge cluster into a light cluster by applying the chargeToPhotonFactor to every point
     *
     *  @param  depositionVector the input charge cluster
     *
     *  @return the output light cluster
     */
    flashana::QCluster_t GetLightCluster(const DepositionVector &depositionVector);

    // Producer labels
    std::string  m_flashLabel;    ///< The label of the flash producer
    std::string  m_pandoraLabel;  ///< The label of the allOutcomes pandora producer

    // Cuts for selecting the beam flash
    float        m_beamWindowStart;  ///< The start time of the beam window
    float        m_beamWindowEnd;    ///< The end time of the beam window
    float        m_minBeamFlashPE;   ///< The minimum number of photoelectrons required to consider a flash as the beam flash

    // Pre-selection cuts to determine if a slice is compatible with the beam flash
    float        m_maxDeltaY;              ///< The maximum difference in Y between the beam flash center and the weighted charge center
    float        m_maxDeltaZ;              ///< The maximum difference in Z between the beam flash center and the weighted charge center
    float        m_maxDeltaYSigma;         ///< As for maxDeltaY, but measured in units of the flash width in Y
    float        m_maxDeltaZSigma;         ///< As for maxDeltaZ, but measured in units of the flash width in Z
    float        m_minChargeToLightRatio;  ///< The minimum ratio between the total charge and the total PE
    float        m_maxChargeToLightRatio;  ///< The maximum ratio between the total charge and the total PE

    // Variables required for flash matching
    float                          m_chargeToNPhotonsTrack;   ///< The conversion factor between charge and number of photons for tracks
    float                          m_chargeToNPhotonsShower;  ///< The conversion factor between charge and number of photons for showers
    flashana::FlashMatchManager    m_flashMatchManager;       ///< The flash match manager
    std::vector<unsigned int>      m_opDetVector;             ///< The ordered vector of optical detector IDs
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
    m_minBeamFlashPE(pset.get<float>("BeamFlashPEThreshold")),
    m_maxDeltaY(pset.get<float>("MaxDeltaY")),
    m_maxDeltaZ(pset.get<float>("MaxDeltaZ")),
    m_maxDeltaYSigma(pset.get<float>("MaxDeltaYSigma")),
    m_maxDeltaZSigma(pset.get<float>("MaxDeltaZSigma")),
    m_minChargeToLightRatio(pset.get<float>("MinChargeToLightRatio")),
    m_maxChargeToLightRatio(pset.get<float>("MaxChargeToLightRatio")),
    m_chargeToNPhotonsTrack(pset.get<float>("ChargeToNPhotonsTrack")),
    m_chargeToNPhotonsShower(pset.get<float>("ChargeToNPhotonsShower"))
{
    m_flashMatchManager.Configure(pset.get<flashana::Config_t>("FlashMatchConfig"));
    
    this->GetOrderedOpDetVector(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::GetOrderedOpDetVector(fhicl::ParameterSet const &pset)
{
    art::ServiceHandle<geo::Geometry> geometry;
    const auto nOpDets(geometry->NOpDets());

    // Get the OpDets in their default order
    std::vector<unsigned int> opDetVector;
    for (unsigned int iChannel = 0; iChannel < nOpDets; ++iChannel)
        opDetVector.push_back(geometry->OpDetFromOpChannel(iChannel));

    // Get the remapped OpDets if required
    if (pset.get<bool>("ShouldRemapPMTs"))
    {
        const auto pmtMapping(pset.get<std::vector<unsigned int> >("OrderedPMTList"));
        
        // Ensure there are the correct number of OpDets
        if (pmtMapping.size() != nOpDets)
            throw cet::exception("FlashNeutrinoId") << "The input PMT remapping vector has the wrong size. Expected " << nOpDets << " elements." << std::endl;

        for (const auto &opDet : pmtMapping)
        {
            // Each OpDet in the default list must be listed once and only once
            if (std::count(opDetVector.begin(), opDetVector.end(), opDet) != 1)
                throw cet::exception("FlashNeutrinoId") << "Unknown or repeated PMT ID: " << opDet << std::endl;

            m_opDetVector.push_back(opDet);
        }
    }
    else
    {
        m_opDetVector = opDetVector;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void FlashNeutrinoId::ClassifySlices(SliceVector &slices, const art::Event &evt) 
{
    if (slices.empty())
        return;

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

    // TODO refactor this into functions
    bool foundViableSlice(false);
    float highestFlashMatchScore(-std::numeric_limits<float>::max());
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        // Collect all spacepoints in the slice that were produced from a hit on the collection plane, and assign them the corresponding charge
        const auto depositionVector(this->GetDepositionVector(pfParticleMap, pfParticleToSpacePointMap, spacePointToHitMap, slices.at(sliceIndex)));

        if (depositionVector.empty())
            continue;

        // Apply the pre-selection cuts to ensure that the slice is compatible with the beam flash
        if (!this->IsSliceCompatibleWithBeamFlash(depositionVector, beamFlash))
            continue;

        // Apply flash-matching, and store this slice if it has the current highest score
        const auto flashMatchScore(this->GetFlashMatchScore(depositionVector, beamFlash));
        if (flashMatchScore > highestFlashMatchScore)
        {
            highestFlashMatchScore = flashMatchScore;
            bestSliceIndex = sliceIndex;
            foundViableSlice = true;
        }
    }

    // Select the best slice (if any)
    if (foundViableSlice)
        slices.at(bestSliceIndex).TagAsTarget();
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

FlashNeutrinoId::DepositionVector FlashNeutrinoId::GetDepositionVector(const PFParticleMap &pfParticleMap,
    const PFParticlesToSpacePoints &pfParticleToSpacePointMap, const SpacePointsToHits &spacePointToHitMap, const Slice &slice) const
{
    // Collect all PFParticles in the slice, including those downstream of the primaries
    // ATTN here we only use the neutrino hypothesis, in theory this should work with either (or indeed both with some thought)
    PFParticleVector allParticlesInSlice;
    this->CollectDownstreamPFParticles(pfParticleMap, slice.GetTargetHypothesis(), allParticlesInSlice);

    DepositionVector depositionVector;
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
            const auto charge(hit->Integral());

            depositionVector.emplace_back(position[0], position[1], position[2], charge, this->GetNPhotons(charge, particle));
        }
    }

    return depositionVector;
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

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &particle) const
{
    return (LArPandoraHelper::IsTrack(particle) ? m_chargeToNPhotonsTrack : m_chargeToNPhotonsShower);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool FlashNeutrinoId::IsSliceCompatibleWithBeamFlash(const FlashNeutrinoId::DepositionVector &depositionVector, const recob::OpFlash &beamFlash) const
{
    // Check the flash is usable
    if (beamFlash.TotalPE() <= std::numeric_limits<float>::epsilon())
        return false;
    
    if (beamFlash.YWidth() <= std::numeric_limits<float>::epsilon())
        return false;
    
    if (beamFlash.ZWidth() <= std::numeric_limits<float>::epsilon())
        return false;

    // Calculate the pre-selection variables
    const auto chargeCenter(this->GetChargeWeightedCenter(depositionVector));
    const auto deltaY(std::abs(chargeCenter.GetY() - beamFlash.YCenter()));
    const auto deltaZ(std::abs(chargeCenter.GetZ() - beamFlash.ZCenter()));
    const auto chargeToLightRatio(this->GetTotalCharge(depositionVector) / beamFlash.TotalPE());  // TODO ATTN check if this should be total PE or max PE. Code differs from technote

    // Check if the slice passes the pre-selection cuts
    return (deltaY < m_maxDeltaY                           &&
            deltaZ < m_maxDeltaZ                           &&
            deltaY / beamFlash.YWidth() < m_maxDeltaYSigma &&
            deltaZ / beamFlash.ZWidth() < m_maxDeltaZSigma &&
            chargeToLightRatio > m_minChargeToLightRatio   &&
            chargeToLightRatio < m_maxChargeToLightRatio   );
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector FlashNeutrinoId::GetChargeWeightedCenter(const FlashNeutrinoId::DepositionVector &depositionVector) const
{
    pandora::CartesianVector center(0.f, 0.f, 0.f);
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
    {
        center += pandora::CartesianVector(chargePoint.m_x, chargePoint.m_y, chargePoint.m_z) * chargePoint.m_charge;
        totalCharge += chargePoint.m_charge;
    }

    if (totalCharge <= std::numeric_limits<float>::epsilon())
        throw cet::exception("FlashNeutrinoId") << "Can't find charge weighted center of slice with zero total charge" << std::endl;

    center *= (1.f / totalCharge);

    return center;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::GetTotalCharge(const FlashNeutrinoId::DepositionVector &depositionVector) const
{
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
        totalCharge += chargePoint.m_charge;

    return totalCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float FlashNeutrinoId::GetFlashMatchScore(const FlashNeutrinoId::DepositionVector &depositionVector, const recob::OpFlash &beamFlash)
{
    m_flashMatchManager.Reset();

    // Convert the flash and the charge cluster into the required format for flash matching
    auto flash(this->ConvertFlashFormat(beamFlash));
    auto lightCluster(this->GetLightCluster(depositionVector));

    // Perform the match
    m_flashMatchManager.Emplace(std::move(flash));
    m_flashMatchManager.Emplace(std::move(lightCluster));
    const auto matches(m_flashMatchManager.Match());

    // Unable to match
    if (matches.empty())
        return 0.f;

    if (matches.size() != 1)
        throw cet::exception("FlashNeutrinoId") << "Flash matching returned multiple matches!" << std::endl;

    return matches.front().score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::Flash_t FlashNeutrinoId::ConvertFlashFormat(const recob::OpFlash &beamFlash) const
{
    // Ensure the input flash is valid
    const auto nOpDets(m_opDetVector.size());
    if (beamFlash.PEs().size() != nOpDets)
        throw cet::exception("FlashNeutrinoId") << "Number of channels in beam flash doesn't match the number of OpDets!" << std::endl;

    // Set the flash properties
    flashana::Flash_t flash;
    flash.x = 0;
    flash.x_err = 0;
    flash.y = beamFlash.YCenter();
    flash.y_err = beamFlash.YWidth();
    flash.z = beamFlash.ZCenter();
    flash.z_err = beamFlash.ZWidth();
    flash.time = beamFlash.Time();
    flash.pe_v.resize(nOpDets);
    flash.pe_err_v.resize(nOpDets);

    // Fill the flash with the PE spectrum
    for (unsigned int i = 0; i < nOpDets; ++i)
    {
        const auto opDet(m_opDetVector.at(i));
        if (opDet < 0 || opDet >= nOpDets)
            throw cet::exception("FlashNeutrinoId") << "OpDet ID, " << opDet << ", is out of range: 0 - " << (nOpDets-1) << std::endl;

        const auto PE(beamFlash.PE(i));
        flash.pe_v.at(opDet) = PE;
        flash.pe_err_v.at(opDet) = std::sqrt(PE);
    }

    return flash;
}

//------------------------------------------------------------------------------------------------------------------------------------------

flashana::QCluster_t FlashNeutrinoId::GetLightCluster(const FlashNeutrinoId::DepositionVector &depositionVector)
{
    flashana::QCluster_t lightCluster;

    for (const auto &point : depositionVector)
        lightCluster.emplace_back(point.m_x, point.m_y, point.m_z, point.m_nPhotons);

    return lightCluster;
}        

//------------------------------------------------------------------------------------------------------------------------------------------

FlashNeutrinoId::Deposition::Deposition(const float x, const float y, const float z, const float charge, const float nPhotons) :
    m_x(x),
    m_y(y),
    m_z(z),
    m_charge(charge),
    m_nPhotons(nPhotons)
{
}

} // namespace lar_pandora
