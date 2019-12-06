/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Tracks
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"

namespace lar_pandora
{

//-----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<art::Ptr<recob::Hit>> LArPandoraTrackUtils::GetHits(const art::Ptr<recob::Track> &track, const art::Event &evt, const std::string &label)
{    
    return LArPandoraTrackUtils::GetAssocProductVector<recob::Hit>(track,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraTrackUtils::GetSpacePoints(const art::Ptr<recob::Track> &track, const art::Event &evt, const std::string &label)
{
    return LArPandoraTrackUtils::GetAssocProductVector<recob::SpacePoint>(track,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const art::Ptr<recob::PFParticle> LArPandoraTrackUtils::GetPFParticle(const art::Ptr<recob::Track> &track, const art::Event &evt, const std::string &label)
{
    return LArPandoraTrackUtils::GetAssocProduct<recob::PFParticle>(track,evt,label,label);
}    

//-----------------------------------------------------------------------------------------------------------------------------------------

const art::Ptr<anab::Calorimetry> LArPandoraTrackUtils::GetCalorimetry(const art::Ptr<recob::Track> &track, const art::Event &evt, const std::string &trackLabel, const std::string &caloLabel)
{
    return LArPandoraTrackUtils::GetAssocProduct<anab::Calorimetry>(track,evt,trackLabel,caloLabel);
}

} // namespace lar_pandora

