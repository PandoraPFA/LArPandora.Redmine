/**
 *  @file   larpandora/LArPandoraAnalysis/AnalysisUtilityExample_module.cc
 *
 *  @brief  This module uses the analysis utilities to demonstrate 
 *          some of their usage. This can be used as a basis for 
 *          writing analysis code using these tools
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraTrackUtils.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraShowerUtils.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  AnalysisUtilityExample class
 */
class AnalysisUtilityExample : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     AnalysisUtilityExample(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~AnalysisUtilityExample();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);

private:

     std::string  m_particleLabel;         ///<
     std::string  m_trackLabel;            ///<
     std::string  m_showerLabel;           ///<
};

DEFINE_ART_MODULE(AnalysisUtilityExample)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/AnalysisBase/T0.h"

#include <iostream>

namespace lar_pandora
{

AnalysisUtilityExample::AnalysisUtilityExample(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset),
    m_particleLabel(pset.get<std::string>("PFParticleModule","pandora")),
    m_trackLabel(pset.get<std::string>("TrackModule","pandoraTrack")),
    m_showerLabel(pset.get<std::string>("ShowerModule","pandoraShower"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AnalysisUtilityExample::~AnalysisUtilityExample()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisUtilityExample::beginJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisUtilityExample::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisUtilityExample::analyze(const art::Event &evt)
{
    // In this example, if we have no neutrino then move on
    const bool hasNeutrino = LArPandoraEventUtils::HasNeutrino(evt,m_particleLabel);
    if (!hasNeutrino)
    {
      return;
    }
    // Now we can get the neutrino
    const art::Ptr<recob::PFParticle> neutrino = LArPandoraEventUtils::GetNeutrino(evt,m_particleLabel);  
    // And its children
    const std::vector<art::Ptr<recob::PFParticle>> neutrinoChildren = LArPandoraPFParticleUtils::GetChildParticles(neutrino,evt,m_particleLabel);
    std::cout << "Found the neutrino, and it has pdg code " << neutrino->PdgCode() << " and " << neutrinoChildren.size() << " child particles" << std::endl;

    // Lets see how many of the children are track- or shower-like
    for (unsigned int iChild = 0; iChild < neutrinoChildren.size(); ++iChild)
    {
        if (LArPandoraPFParticleUtils::IsTrack(neutrinoChildren.at(iChild),evt,m_particleLabel,m_trackLabel))
        {
            std::cout << "Child " << iChild << " is track-like" << std::endl;
        }
        else if (LArPandoraPFParticleUtils::IsShower(neutrinoChildren.at(iChild),evt,m_particleLabel,m_showerLabel))
        {
            std::cout << "Child " << iChild << " is shower-like" << std::endl;
        }
        else
        {
            std::cout << "Child " << iChild << " has no track or shower association" << std::endl;   
        }
    }

    // Look for the vertex
    const art::Ptr<recob::Vertex> neutrinoVertex = LArPandoraPFParticleUtils::GetVertex(neutrino,evt,m_particleLabel);
    std::cout << "Interaction vertex at " << neutrinoVertex->position().X() << ", "
                                          << neutrinoVertex->position().Y() << ", "
                                          << neutrinoVertex->position().Z() << std::endl;
  
    // Let's have a look at all of the particles
    const std::vector<art::Ptr<recob::PFParticle>> recoParticles = LArPandoraEventUtils::GetPFParticles(evt,m_particleLabel);

    unsigned int nTracks = 0;
    unsigned int nShowers = 0;
    unsigned int nHits = 0;
    for (const art::Ptr<recob::PFParticle> &p : recoParticles)
    {
        if (LArPandoraPFParticleUtils::IsTrack(p,evt,m_particleLabel,m_trackLabel))
        {
            const art::Ptr<recob::Track> track = LArPandoraPFParticleUtils::GetTrack(p,evt,m_particleLabel,m_trackLabel);
            nHits += LArPandoraTrackUtils::GetHits(track,evt,m_trackLabel).size();
            ++nTracks;
        }
        else if (LArPandoraPFParticleUtils::IsShower(p,evt,m_particleLabel,m_showerLabel))
        {
            const art::Ptr<recob::Shower> shower = LArPandoraPFParticleUtils::GetShower(p,evt,m_particleLabel,m_showerLabel);
            nHits += LArPandoraShowerUtils::GetHits(shower,evt,m_showerLabel).size();
            ++nShowers;
        }
        else
        {
            nHits += LArPandoraPFParticleUtils::GetHits(p,evt,m_particleLabel).size();
        }
    }
    std::cout << "Total number of hits in all " << recoParticles.size() << " PFParticles = " << nHits << std::endl; 
}

} //namespace lar_pandora

