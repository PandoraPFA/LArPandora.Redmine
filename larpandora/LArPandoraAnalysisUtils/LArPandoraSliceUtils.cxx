/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraSliceUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Slices
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraSliceUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

namespace lar_pandora
{

const std::vector<art::Ptr<recob::Hit>> LArPandoraSliceUtils::GetHits(const art::Ptr<recob::Slice> &slice, const art::Event &evt, const std::string &label)
{    
    return GetAssocProductVector<recob::Hit>(slice,evt,label,label);
}

} // namespace lar_pandora


