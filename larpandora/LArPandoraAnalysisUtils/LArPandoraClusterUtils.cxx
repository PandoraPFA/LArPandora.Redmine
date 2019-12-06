/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraClusterUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Clusters
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraClusterUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace lar_pandora
{

const std::vector<art::Ptr<recob::Hit>> LArPandoraClusterUtils::GetHits(const art::Ptr<recob::Cluster> &cluster, const art::Event &evt, const std::string &label)
{    
    return GetAssocProductVector<recob::Hit>(cluster,evt,label,label);
}


} // namespace lar_pandora


