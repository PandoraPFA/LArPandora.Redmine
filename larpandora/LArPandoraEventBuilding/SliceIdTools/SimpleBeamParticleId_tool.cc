/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraSimpleBeamParticleId_tool.cc
 *
 *  @brief  implementation of the lar pandora simple beam particle id tool
 */

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "larpandora/LArPandoraEventBuilding/SliceIdBaseTool.h"
#include "larpandora/LArPandoraEventBuilding/Slice.h"

namespace lar_pandora
{

/**
 *  @brief  Simple beam particle ID tool that selects the most likely beam particle slice using the scores from Pandora
 */
class SimpleBeamParticleId : SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  pset FHiCL parameter set
     */
    SimpleBeamParticleId(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Classify slices as neutrino or cosmic
     *
     *  @param  slices the input vector of slices to classify
     *  @param  evt the art event
     */
    void ClassifySlices(SliceVector &slices, const art::Event &evt) override;
    
private:
    float m_minBDTScore; ///< The minimum BDT score to select a slice as a beam particle

};

DEFINE_ART_CLASS_TOOL(SimpleBeamParticleId)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace lar_pandora
{
 
SimpleBeamParticleId::SimpleBeamParticleId(fhicl::ParameterSet const &pset) :
    m_minBDTScore(pset.get<float>("MinBDTScore"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleBeamParticleId::ClassifySlices(SliceVector &slices, const art::Event &/*evt*/) 
{
    if (slices.empty()) return;

    // Find the most probable slice
    float highestBeamParticleScore(-std::numeric_limits<float>::max());
    unsigned int mostProbableSliceIndex(std::numeric_limits<unsigned int>::max());

    for (unsigned int sliceIndex = 0; sliceIndex < slices.size(); ++sliceIndex)
    {
        const float beamParticleScore(slices.at(sliceIndex).GetTopologicalScore());
        std::cout << "Slice " << sliceIndex << " - " << beamParticleScore << std::endl;
        if (beamParticleScore > highestBeamParticleScore)
        {
            highestBeamParticleScore = beamParticleScore;
            mostProbableSliceIndex = sliceIndex;
        }
    }

    std::cout << "Tagging slice " << mostProbableSliceIndex << std::endl;

    // Tag the most probable slice as a the beam particle
    slices.at(mostProbableSliceIndex).TagAsTarget();
}

} // namespace lar_pandora
