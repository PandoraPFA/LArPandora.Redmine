/**
 *  @file   larpandora/LArPandoraEventBuilding/Slice.h
 *
 *  @brief  header for the lar pandora slice class
 */

#ifndef LAR_PANDORA_NEUTRINO_ID_BASE_TOOL_H
#define LAR_PANDORA_NEUTRINO_ID_BASE_TOOL_H 1

#include "larpandora/LArPandoraEventBuilding/Slice.h"

namespace lar_pandora
{

/**
 *  @brief  Abstract base class for a neutrino ID tool
 */
class NeutrinoIdBaseTool
{
public:
    virtual ~NeutrinoIdBaseTool() noexcept = default;

    /**
     *  @brief  The tools interface function. Here the derived tool will classify the input slices
     *
     *  @param  slices the input vector of slices to classify
     */
    virtual void ClassifySlices(SliceVector &slices) = 0;
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_NEUTRINO_ID_BASE_TOOL_H
