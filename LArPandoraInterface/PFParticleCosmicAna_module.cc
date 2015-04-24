/**
 *  @file   larpandora/LArPandoraInterface/PFParticleCosmicAna_module.cc
 *
 *  @brief  Analysis module for created particles
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

// ROOT includes
#include "TTree.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  PFParticleCosmicAna class
 */
class PFParticleCosmicAna : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
     PFParticleCosmicAna(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~PFParticleCosmicAna();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);
     void reconfigure(fhicl::ParameterSet const &pset);

private:

     TTree       *m_pCosmicTree;         ///< 

     int          m_run;                   ///< 
     int          m_event;                 ///< 
     int          m_index;                 ///<

     int          m_self;                  ///<
     int          m_pdgCode;               ///<
     int          m_isTrack;               ///<
     int          m_isPrimary;             ///<
     float        m_cosmicScore;           ///<
     int          m_cosmicTag;             ///<
     int          m_nTracks;               ///<
     int          m_nHits;                 ///<
     int          m_nCosmicHits;           ///<
     int          m_nNeutrinoHits;         ///<     

     std::string  m_hitfinderLabel;         ///<
     std::string  m_trackfitLabel;          ///<
     std::string  m_particleLabel;          ///<
     std::string  m_cosmicLabel;            ///<
     std::string  m_geantModuleLabel;       ///<

     bool         m_useDaughterPFParticles; ///<
     bool         m_useDaughterMCParticles; ///<

     float        m_cosmicScoreThreshold;   ///<
};

DEFINE_ART_MODULE(PFParticleCosmicAna)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
/**
 *  @file   LArPandora/PFParticleCosmicAna.cxx
 * 
 *  @brief  Implementation of PFParticle analysis module
 * 
 *  $Log: $
 */

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FindOneP.h"

// LArSoft includes
#include "AnalysisBase/CosmicTag.h"

// Local includes
#include "LArPandoraInterface/LArPandoraCollector.h"

// std includes
#include <iostream>

namespace lar_pandora
{

PFParticleCosmicAna::PFParticleCosmicAna(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleCosmicAna::~PFParticleCosmicAna()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleCosmicAna::reconfigure(fhicl::ParameterSet const &pset)
{ 
    m_cosmicLabel = pset.get<std::string>("CosmicTagModule","cosmictagger");
    m_particleLabel = pset.get<std::string>("PFParticleModule","pandora");
    m_trackfitLabel = pset.get<std::string>("TrackFitModule","trackfit");
    m_hitfinderLabel = pset.get<std::string>("HitFinderModule","gaushit");
    m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");

    m_useDaughterPFParticles = pset.get<bool>("UseDaughterPFParticles",true);
    m_useDaughterMCParticles = pset.get<bool>("UseDaughterMCParticles",true);

    m_cosmicScoreThreshold = pset.get<float>("CosmicScoreThreshold",0.75);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleCosmicAna::beginJob()
{
    mf::LogDebug("LArPandora") << " *** PFParticleCosmicAna::beginJob() *** " << std::endl; 

    // 
    art::ServiceHandle<art::TFileService> tfs;
 
    m_pCosmicTree = tfs->make<TTree>("pandora", "LAr Cosmic Tree");
    m_pCosmicTree->Branch("run", &m_run, "run/I");
    m_pCosmicTree->Branch("event", &m_event, "event/I");
    m_pCosmicTree->Branch("index", &m_index, "index/I");
    m_pCosmicTree->Branch("self", &m_self, "self/I");
    m_pCosmicTree->Branch("pdgCode", &m_pdgCode, "pdgCode/I"); 
    m_pCosmicTree->Branch("isTrack", &m_isTrack, "isTrack/I");
    m_pCosmicTree->Branch("isPrimary", &m_isPrimary, "isPrimary/I");
    m_pCosmicTree->Branch("cosmicTag", &m_cosmicTag, "cosmicTag/I");
    m_pCosmicTree->Branch("cosmicScore", &m_cosmicScore, "cosmicScore/F");
    m_pCosmicTree->Branch("nTracks", &m_nTracks, "nTracks/I");
    m_pCosmicTree->Branch("nHits", &m_nHits, "nHits/I");
    m_pCosmicTree->Branch("nCosmicHits", &m_nCosmicHits, "nCosmicHits/I");
    m_pCosmicTree->Branch("nNeutrinoHits", &m_nNeutrinoHits, "nNeutrinoHits/I");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleCosmicAna::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleCosmicAna::analyze(const art::Event &evt)
{
    std::cout << " *** PFParticleCosmicAna::analyze(...) *** " << std::endl;

    m_run = evt.run();
    m_event = evt.id().event();
    m_index = 0;

    m_self = 0;
    m_pdgCode = 0;
    m_isTrack = 0;
    m_isPrimary = 0;
    m_cosmicTag = 0;
    m_cosmicScore = 0.f;
    m_nTracks = 0;
    m_nHits = 0;
    m_nCosmicHits = 0;
    m_nNeutrinoHits = 0;

    std::cout << "  Run: " << m_run << std::endl;
    std::cout << "  Event: " << m_event << std::endl; 

    // Collect True Particles
    // ======================
    HitVector hitVector;
    MCTruthToMCParticles truthToParticles;
    MCParticlesToMCTruth particlesToTruth;
    MCParticlesToHits trueParticlesToHits;
    HitsToMCParticles trueHitsToParticles;

    LArPandoraCollector::CollectHits(evt, m_hitfinderLabel, hitVector);
    LArPandoraCollector::CollectMCParticles(evt, m_geantModuleLabel, truthToParticles, particlesToTruth);
    LArPandoraCollector::BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles,
        (m_useDaughterMCParticles ? LArPandoraCollector::kAddDaughters : LArPandoraCollector::kIgnoreDaughters));


    // Collect Reco Particles
    // ======================
    PFParticleVector recoParticleVector;
    PFParticlesToHits recoParticlesToHits;
    HitsToPFParticles recoHitsToParticles;

    LArPandoraCollector::CollectPFParticles(evt, m_particleLabel, recoParticleVector);
    LArPandoraCollector::BuildPFParticleHitMaps(evt, m_particleLabel, m_particleLabel, recoParticlesToHits, recoHitsToParticles, 
        (m_useDaughterPFParticles ? LArPandoraCollector::kAddDaughters : LArPandoraCollector::kIgnoreDaughters));

    std::cout << "  PFParticles: " << recoParticleVector.size() << std::endl;


    // Collect Reco Tracks
    // ===================
    TrackVector recoTrackVector;
    PFParticlesToTracks recoParticlesToTracks;
    LArPandoraCollector::CollectTracks(evt, m_trackfitLabel, recoTrackVector, recoParticlesToTracks);


    // Collect Cosmic Tags
    // ===================
    typedef std::map<art::Ptr<recob::Track>, art::Ptr<anab::CosmicTag>> TracksToCosmicTags;
    TracksToCosmicTags recoTracksToCosmicTags;

    art::Handle< std::vector<anab::CosmicTag> > theCosmicTags;
    evt.getByLabel(m_cosmicLabel, theCosmicTags); // Note: in general, there could be many tagging algorithms

    if (theCosmicTags.isValid())
    {
        art::FindOneP<recob::Track> theCosmicAssns(theCosmicTags, evt, m_cosmicLabel);
        for (unsigned int i = 0; i < theCosmicTags->size(); ++i)
        {
            const art::Ptr<anab::CosmicTag> cosmicTag(theCosmicTags, i);
            const art::Ptr<recob::Track> track = theCosmicAssns.at(i);
            recoTracksToCosmicTags[track] = cosmicTag;
        }
    }


    // Loop over Reco Particles
    // ========================
    for (PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        const HitVector &hitVector = iter1->second;

        if (hitVector.empty())
            continue;
 
        m_self          = recoParticle->Self();
        m_pdgCode       = recoParticle->PdgCode();
        m_isPrimary     = recoParticle->IsPrimary();
        m_isTrack       = LArPandoraCollector::IsTrack(recoParticle);
        m_cosmicScore   = 0.f;
        m_cosmicTag     = 0; 
        m_nTracks       = 0; 
        m_nHits         = hitVector.size();
        m_nNeutrinoHits = 0;
        m_nCosmicHits   = 0;

        // Get cosmic tags associated with this particle
        PFParticlesToTracks::const_iterator iter2 = recoParticlesToTracks.find(recoParticle);
        if (recoParticlesToTracks.end() != iter2)
        {
            for (TrackVector::const_iterator iter3 = iter2->second.begin(), iterEnd3 = iter2->second.end(); iter3 != iterEnd3; ++iter3)
            {
                const art::Ptr<recob::Track> track = *iter3;

                TracksToCosmicTags::const_iterator iter4 = recoTracksToCosmicTags.find(track);
                if (recoTracksToCosmicTags.end() != iter4)
                {
                    const art::Ptr<anab::CosmicTag> cosmicTag = iter4->second;
                    const float thisScore(cosmicTag->CosmicScore());
                    if (thisScore > m_cosmicScore)
                        m_cosmicScore = thisScore;
                }

                ++m_nTracks;
            }

            if (m_cosmicScore > m_cosmicScoreThreshold)
                m_cosmicTag = 1; 
        }

        // Loop over hits associated with this particle
        for (HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

            MCParticlesToMCTruth::const_iterator iter4 = particlesToTruth.find(trueParticle);
            if (particlesToTruth.end() == iter4)
                throw cet::exception("LArPandora") << " PFParticleCosmicAna::analyze --- Found a true particle without any ancestry information ";
        
            const art::Ptr<simb::MCTruth> truth = iter4->second;

            if (truth->NeutrinoSet())
            {
                ++m_nNeutrinoHits;
            }
                else
            {
                ++m_nCosmicHits;
            }
        }

        std::cout << "   PFParticle: [" << m_index << "] nCosmicHits=" << m_nCosmicHits << ", nNeutrinoHits=" << m_nNeutrinoHits << ", nTracks=" << m_nTracks << ", cosmicTag=" << m_cosmicTag << ", cosmicScore=" << m_cosmicScore << std::endl;

        m_pCosmicTree->Fill();
        ++m_index;
    }
}

} //namespace lar_pandora
