/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdTools/FlashNeutrinoId_tool.cc
 *
 *  @brief  implementation of the flash based neutrino id tool
 */

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/OpFlash.h"

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

    std::string  m_flashLabel;       ///< The label of the flash producer
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
    m_beamWindowStart(pset.get<float>("BeamWindowStartTime", 3.2f)),
    m_beamWindowEnd(pset.get<float>("BeamWindowEndTime", 4.8f)),
    m_minBeamFlashPE(pset.get<float>("BeamFlashPEThreshold", 50.f))
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
   
    /*
    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
    }*/
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

} // namespace lar_pandora
