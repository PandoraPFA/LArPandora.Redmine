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

    const std::vector<art::Ptr<anab::T0>> LArPandoraPFParticleUtils::GetT0(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<anab::T0>> theseT0s;
        GetAssocProductVector(part,evt,label,label,theseT0s);
        return theseT0s;
    }

    const std::vector<art::Ptr<anab::CosmicTag>> LArPandoraPFParticleUtils::GetCosmicTag(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<anab::CosmicTag>> theseTags;
        GetAssocProductVector(part,evt,label,label,theseTags); 
        return theseTags;
    }

    const std::vector<art::Ptr<recob::PFParticle>> LArPandoraPFParticleUtils::GetChildParticles(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {

        art::Handle<std::vector<recob::PFParticle>> particles;
        evt.getByLabel(label,particles);

        unsigned int thisParticleIndex = part.key();

        std::vector<art::Ptr<recob::PFParticle>> children;

        for (unsigned int p = 0; p < particles->size(); ++p)
        {     
            if (particles->at(p).Parent() == thisParticleIndex)
            {
                art::Ptr<recob::PFParticle> pPtr(particles,p);
                children.push_back(pPtr);
            }
        }

        return children;
    }

    const std::vector<art::Ptr<recob::Hit>> LArPandoraPFParticleUtils::GetHits(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {    
        // There isn't a direct association between PFParticles and hits, so we go via clusters
        std::vector<art::Ptr<recob::Cluster>> theseClusters;
        GetAssocProductVector(part,evt,label,label,theseClusters);

        std::vector<art::Ptr<recob::Hit>> theseHits;
        for (const art::Ptr<recob::Cluster> cluster : theseClusters)
        {
          std::vector<art::Ptr<recob::Hit>> tempHits;
          GetAssocProductVector(cluster,evt,label,label,tempHits);
          theseHits.insert(theseHits.end(),tempHits.begin(),tempHits.end());
        }
        return theseHits;
    }

    const std::vector<art::Ptr<recob::SpacePoint>> LArPandoraPFParticleUtils::GetSpacePoints(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {

        std::vector<art::Ptr<recob::SpacePoint>> theseSPs;
        GetAssocProductVector(part,evt,label,label,theseSPs);
        return theseSPs;
    }

    const art::Ptr<recob::Track> LArPandoraPFParticleUtils::GetTrack(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &particleLabel, const std::string &trackLabel)
    {
        std::vector<art::Ptr<recob::Track>> theseTracks;
        GetAssocProductVector(part,evt,particleLabel,trackLabel,theseTracks);

        if (theseTracks.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraPFParticleUtils::GetTrack --- No associated track found";
        }
        return theseTracks.at(0);
    }    
    
    const art::Ptr<recob::Shower> LArPandoraPFParticleUtils::GetShower(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &particleLabel, const std::string &showerLabel)
    {

        std::vector<art::Ptr<recob::Shower>> theseShowers;
        GetAssocProductVector(part,evt,particleLabel,showerLabel,theseShowers);

        if (theseShowers.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraPFParticleUtils::GetShower --- No associated shower found";
        }
        return theseShowers.at(0);

    }

    const art::Ptr<recob::Vertex> LArPandoraPFParticleUtils::GetVertex(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &particleLabel)
    {

        std::vector<art::Ptr<recob::Vertex>> theseVertices;
        GetAssocProductVector(part,evt,particleLabel,particleLabel,theseVertices);

        if (theseVertices.size() == 0)
        {
            throw cet::exception("LArPandora") << "LArPandoraPFParticleUtils::GetVertex --- No associated vertex found";
        }
        return theseVertices.at(0);
    }

    const art::Ptr<recob::Slice> LArPandoraPFParticleUtils::GetSlice(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {
        const art::Ptr<larpandoraobj::PFParticleMetadata> thisMetadata = GetMetadata(part,evt,label);
        std::map<std::string,float> metaMap = thisMetadata->GetPropertiesMap();

        unsigned int sliceIndex;

        if (metaMap.find("SliceIndex") == metaMap.end())   
        {
            throw cet::exception("LArPandora") << " LArPandoraPFParticleUtils::GetSlice --- No associated slice found";
        }
        else{
            sliceIndex = static_cast<unsigned int>(metaMap.at("SliceIndex"));
        }

        std::vector<art::Ptr<recob::Slice>> theseSlices;
        GetProductVector(evt,label,theseSlices);

        return theseSlices.at(sliceIndex);
    }

    const art::Ptr<larpandoraobj::PFParticleMetadata> LArPandoraPFParticleUtils::GetMetadata(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &label)
    {
        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> theseMetadata;
        GetAssocProductVector(part,evt,label,label,theseMetadata);

        if (theseMetadata.size() == 0)
        {
            throw cet::exception("LArPandora") << " LArPandoraPFParticleUtils::GetMetadata --- No associated metadata found";
        }

        return theseMetadata.at(0);
    }

    bool LArPandoraPFParticleUtils::IsTrack(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &particleLabel, const std::string &trackLabel)
    {
        // This function needs to fail if GetTrack would fail
        std::vector<art::Ptr<recob::Track>> theseTracks;
        GetAssocProductVector(part,evt,particleLabel,trackLabel,theseTracks);
        if (theseTracks.size() == 0)
        {
            return false;
        }
        else return true;
    }

    bool LArPandoraPFParticleUtils::IsShower(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &particleLabel, const std::string &showerLabel)
    {
        std::vector<art::Ptr<recob::Shower>> theseShowers;
        GetAssocProductVector(part,evt,particleLabel,showerLabel,theseShowers);

        if (theseShowers.size() == 0)
        {
            return false;
        }
        else return true;
    }

    bool LArPandoraPFParticleUtils::IsClearCosmic(const art::Ptr<recob::PFParticle> &part, art::Event const &evt, const std::string &particleLabel)
    {
        const art::Ptr<larpandoraobj::PFParticleMetadata> metadata = GetMetadata(part,evt,particleLabel);

        std::map<std::string,float> metaMap = metadata->GetPropertiesMap();

        if (metaMap.find("IsClearCosmic") != metaMap.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    bool LArPandoraPFParticleUtils::IsNeutrino(const art::Ptr<recob::PFParticle> &particle)
    {
        return LArPandoraHelper::IsNeutrino(particle);
    }


} // namespace lar_pandora


