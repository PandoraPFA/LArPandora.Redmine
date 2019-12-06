/**
*
* @file larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"

namespace lar_pandora
{

std::vector<art::Ptr<anab::T0>> LArPandoraPFParticleUtils::GetT0(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return LArPandoraPFParticleUtils::GetAssocProductVector<anab::T0>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<anab::CosmicTag>> LArPandoraPFParticleUtils::GetCosmicTag(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return LArPandoraPFParticleUtils::GetAssocProductVector<anab::CosmicTag>(pParticle,evt,label,label); 
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::PFParticle>> LArPandoraPFParticleUtils::GetChildParticles(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    art::Handle<std::vector<recob::PFParticle>> pParticles;
    bool success = evt.getByLabel(label,pParticles);
    
    if (!success)
    {   
        mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<recob::PFParticle>>();
    }

    std::vector<art::Ptr<recob::PFParticle>> children;

    for (unsigned int index = 0; index < pParticles->size(); ++index)
    {     
        if (pParticles->at(index).Parent() == pParticle.key())
        {
            art::Ptr<recob::PFParticle> pChild(pParticles,index);
            children.push_back(pChild);
        }
    }

    return children;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> LArPandoraPFParticleUtils::GetHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{    
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<art::Ptr<recob::Cluster>> theseClusters = LArPandoraPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<art::Ptr<recob::Hit>> theseHits;
    for (const art::Ptr<recob::Cluster> pCluster : theseClusters)
    {
      const std::vector<art::Ptr<recob::Hit>> tempHits = LArPandoraPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
      theseHits.insert(theseHits.end(),tempHits.begin(),tempHits.end());
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> LArPandoraPFParticleUtils::GetViewHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label, const unsigned short &view)
{
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<art::Ptr<recob::Cluster>> theseClusters = LArPandoraPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<art::Ptr<recob::Hit>> theseHits;
    for (const art::Ptr<recob::Cluster> pCluster : theseClusters)
    { 
        const std::vector<art::Ptr<recob::Hit>> tempHits = LArPandoraPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
        for(const art::Ptr<recob::Hit> pHit : tempHits)
        {
            if (pHit->View() == view)
            {
                theseHits.push_back(pHit);
            }
        }
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> LArPandoraPFParticleUtils::GetSpacePoints(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return LArPandoraPFParticleUtils::GetAssocProductVector<recob::SpacePoint>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Track> LArPandoraPFParticleUtils::GetTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    return LArPandoraPFParticleUtils::GetAssocProduct<recob::Track>(pParticle,evt,particleLabel,trackLabel);
}    

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Shower> LArPandoraPFParticleUtils::GetShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    return LArPandoraPFParticleUtils::GetAssocProduct<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Vertex> LArPandoraPFParticleUtils::GetVertex(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    return LArPandoraPFParticleUtils::GetAssocProduct<recob::Vertex>(pParticle,evt,particleLabel,particleLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Slice> LArPandoraPFParticleUtils::GetSlice(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> pMetadata = LArPandoraPFParticleUtils::GetMetadata(pParticle,evt,label);
    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    unsigned int sliceIndex;
    std::map<std::string,float>::iterator mapItr = metaMap.find("SliceIndex");

    if (mapItr == metaMap.end())   
    {
        throw cet::exception("LArPandora") << " LArPandoraPFParticleUtils::GetSlice --- No associated slice found";
    }
    else
    {
        sliceIndex = mapItr->second;
    }

    return LArPandoraPFParticleUtils::GetProductVector<recob::Slice>(evt,label).at(sliceIndex);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<larpandoraobj::PFParticleMetadata> LArPandoraPFParticleUtils::GetMetadata(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return LArPandoraPFParticleUtils::GetAssocProduct<larpandoraobj::PFParticleMetadata>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraPFParticleUtils::IsTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    // This function needs to fail if GetTrack would fail
    const std::vector<art::Ptr<recob::Track>> theseTracks = LArPandoraPFParticleUtils::GetAssocProductVector<recob::Track>(pParticle,evt,particleLabel,trackLabel);

    return !theseTracks.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraPFParticleUtils::IsShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    const std::vector<art::Ptr<recob::Shower>> theseShowers = LArPandoraPFParticleUtils::GetAssocProductVector<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
   
    return !theseShowers.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraPFParticleUtils::IsClearCosmic(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> pMetadata = LArPandoraPFParticleUtils::GetMetadata(pParticle,evt,particleLabel);

    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    return metaMap.find("IsClearCosmic") != metaMap.end();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArPandoraPFParticleUtils::IsNeutrino(const art::Ptr<recob::PFParticle> &pParticle)
{
    return LArPandoraHelper::IsNeutrino(pParticle);
}

} // namespace lar_pandora


