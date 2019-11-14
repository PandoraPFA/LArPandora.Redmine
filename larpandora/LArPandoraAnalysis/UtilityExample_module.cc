/**
 *  @file   larpandora/LArPandoraAnalysis/UtilityExample_module.cc
 *
 *  @brief  Analysis module for created particles
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "larpandora/LArPandoraAnalysisUtils/LArPandoraEventUtils.h"
#include "larpandora/LArPandoraAnalysisUtils/LArPandoraPFParticleUtils.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  UtilityExample class
 */
class UtilityExample : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     UtilityExample(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~UtilityExample();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);
     void reconfigure(fhicl::ParameterSet const &pset);

private:

     std::string  m_particleLabel;         ///<
     std::string  m_trackLabel;            ///<
     std::string  m_showerLabel;           ///<
};

DEFINE_ART_MODULE(UtilityExample)

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

UtilityExample::UtilityExample(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

UtilityExample::~UtilityExample()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UtilityExample::reconfigure(fhicl::ParameterSet const &pset)
{
    m_particleLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_trackLabel = pset.get<std::string>("TrackModule","pandoraTrack");
    m_showerLabel = pset.get<std::string>("ShowerModule","pandoraShower");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UtilityExample::beginJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UtilityExample::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UtilityExample::analyze(const art::Event &evt)
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
  
    // Start by getting the PFParticles
 //   const std::vector<art::Ptr<recob::PFParticle>> recoParticles = LArPandoraEventUtils::GetPFParticles(evt,m_particleLabel);


}

} //namespace lar_pandora

