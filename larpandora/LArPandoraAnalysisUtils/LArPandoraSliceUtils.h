/**
 *
 * @file larpandora/LArPandoraAnalysisUtils/LArPandoraSliceUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Slices
*/

#ifndef LAR_PANDORA_SLICE_UTILS_H
#define LAR_PANDORA_SLICE_UTILS_H

#include "art/Framework/Principal/Event.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Access the type defs defined in the helper
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraUtilsBase.h"

#include <string>
#include <vector>

namespace lar_pandora
{
/**
 *
 * @brief LArPandoraSliceUtils class
 *
*/
class LArPandoraSliceUtils:LArPandoraUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the slice.
    *
    * @param slice is the slice for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the slice producer
    * 
    * @return vector of art::Ptrs to the hits
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Slice> &pSlice, const art::Event &evt, const std::string &label);
};

} // namespace lar_pandora

#endif // LAR_PANDORA_SLICE_UTILS_H
