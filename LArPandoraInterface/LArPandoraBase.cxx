// Framework includes
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes 
#include "Geometry/Geometry.h"
#include "SimulationBase/MCTruth.h"

// Pandora includes
#include "Objects/ParticleFlowObject.h"

// Local includes
#include "LArContent.h"
#include "LArPandoraBase.h"

// System includes
#include <iostream>

namespace lar_pandora {

LArPandoraBase::LArPandoraBase(fhicl::ParameterSet const &pset)
{
    this->reconfigure(pset);
    m_pPandora = new pandora::Pandora();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraBase::~LArPandoraBase()
{
    delete m_pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::reconfigure(fhicl::ParameterSet const &pset)
{
    m_enableProduction = pset.get<bool>("EnableProduction",true);
    m_enableMCParticles = pset.get<bool>("EnableMCParticles",false);
    m_enableMonitoring = pset.get<bool>("EnableMonitoring",false);

    m_configFile = pset.get<std::string>("ConfigFile");
    m_geantModuleLabel = pset.get<std::string>("GeantModuleLabel","largeant");
    m_hitfinderModuleLabel = pset.get<std::string>("HitFinderModuleLabel","ffthit");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::beginJob()
{
    this->InitializePandora();  

    if (m_enableMonitoring)
        this->InitializeMonitoring();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::produce(art::Event &evt)
{ 
    mf::LogInfo("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] *** " << std::endl;

    cet::cpu_timer theClock;

    HitVector theArtHits;
    HitToParticleMap theHitToParticleMap;
    TruthToParticleMap theTruthToParticleMap;
    ParticleMap theParticleMap;
    HitMap thePandoraHits;

    this->PrepareEvent(evt);
    this->CollectArtHits(evt, theArtHits, theHitToParticleMap);
    this->CreatePandoraHits(theArtHits, thePandoraHits);

    if (m_enableMCParticles && !evt.isRealData())
    {
        this->CollectArtParticles(evt, theParticleMap, theTruthToParticleMap);
        this->CreatePandoraParticles(theParticleMap, theTruthToParticleMap);
        this->CreatePandoraLinks(thePandoraHits, theHitToParticleMap);
    }

    if (m_enableMonitoring)
        theClock.start();

    this->RunPandora();

    if (m_enableMonitoring)
        theClock.stop();

    if (m_enableProduction)
        this->ProduceArtOutput(evt, thePandoraHits);

    this->ResetPandora();

    if (m_enableMonitoring)
    {
        m_run   = evt.run();
        m_event = evt.id().event();
        m_time  = theClock.accumulated_real_time();
        m_hits  = static_cast<int>(theArtHits.size());
        m_pRecoTree->Fill();
    }
   
    mf::LogDebug("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "]  Done! *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::InitializePandora()
{ 
    mf::LogDebug("LArPandora") << " *** LArPandora::InitializePandora(...) *** " << std::endl;

    // Find the Pandora settings file (must be within 'FW_SEARCH_PATH')
    cet::search_path sp("FW_SEARCH_PATH");
    std::string configFileName("");

    mf::LogDebug("LArPandora") << "   Load Pandora settings: " << m_configFile << std::endl;
    mf::LogDebug("LArPandora") << "   Search path: " << sp.to_string() << std::endl;

    if (!sp.find_file(m_configFile, configFileName))
    {
        mf::LogError("LArPandora") << "   Failed to find: " << m_configFile << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }
    else
    {
        mf::LogDebug("LArPandora") << "   Found it: " <<  configFileName << std::endl;
    }
    
    this->CreatePandoraGeometry();

    // Register the algorithms and read the settings
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*m_pPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterHelperFunctions(*m_pPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterResetFunctions(*m_pPandora));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, configFileName));

    mf::LogDebug("LArPandora") << " *** LArPandoraBase::InitializePandora(...)  Done! *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("hits", &m_hits, "hits/I");
    m_pRecoTree->Branch("time", &m_time, "time/F");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::PrepareEvent(const art::Event &evt)
{
    m_run   = evt.run();
    m_event = evt.id().event();
    m_hits  = 0;
    m_time  = 0.f;
 
    if (m_enableMCParticles && !evt.isRealData())
    {
        art::ServiceHandle<cheat::BackTracker> theBackTracker; 

	// Bail out if there is no back-tracking information
        if( theBackTracker->GetSetOfTrackIDs().size() == 0 )
	{
	    mf::LogError("LArPandora") << "   Failed to load back-tracking data " << std::endl;
	    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);  
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::CollectArtHits(const art::Event &evt, HitVector &hitVector, HitToParticleMap &hitToParticleMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandoraBase::CollectArtHits(...) *** " << std::endl;

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(m_hitfinderModuleLabel, hitHandle);

    art::ServiceHandle<cheat::BackTracker> theBackTracker; 

    for (unsigned int iHit = 0, iHitEnd = hitHandle->size(); iHit < iHitEnd; ++iHit)
    {
        art::Ptr<recob::Hit> hit(hitHandle, iHit);
        hitVector.push_back(hit);

        if (m_enableMCParticles && !evt.isRealData())
        {
            const std::vector<cheat::TrackIDE> trackCollection(theBackTracker->HitToTrackID(hit));

            for (unsigned int iTrack = 0, iTrackEnd = trackCollection.size(); iTrack < iTrackEnd; ++iTrack)
            {
                cheat::TrackIDE trackIDE = trackCollection.at(iTrack);
                hitToParticleMap[hit].push_back(trackIDE);
	    }
	}
    }

    mf::LogDebug("LArPandora") << "   Number of ART hits: " << hitVector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::CollectArtParticles(const art::Event &evt, ParticleMap &particleMap, TruthToParticleMap &truthToParticleMap) const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::CollectArtParticles(...) *** " << std::endl;

    art::Handle< std::vector<simb::MCParticle> > mcParticleHandle;
    evt.getByLabel(m_geantModuleLabel, mcParticleHandle);

    art::ServiceHandle<cheat::BackTracker> theBackTracker; 

    for (unsigned int i = 0, iEnd = mcParticleHandle->size(); i < iEnd; ++i)
    {
        art::Ptr<simb::MCParticle> particle(mcParticleHandle, i);
        particleMap[particle->TrackId()] = particle;

        art::Ptr<simb::MCTruth> truth(theBackTracker->TrackIDToMCTruth(particle->TrackId()));
        truthToParticleMap[truth].push_back(particle->TrackId());
    }

    mf::LogDebug("LArPandora") << "   Number of ART particles: " << particleMap.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::RunPandora() const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::RunPandora() *** " << std::endl;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::ResetPandora() const
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ResetPandora() *** " << std::endl;

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraBase::GetStartAndEndPoints(const art::Ptr<simb::MCParticle> &particle, int& startT, int& endT) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    bool foundStartPosition(false);

    const int numTrajectoryPoints(static_cast<int>(particle->NumberTrajectoryPoints()));

    for (int nt = 0; nt < numTrajectoryPoints; ++nt)
    {
        try{
            double pos[3] = {particle->Vx(nt), particle->Vy(nt), particle->Vz(nt)};
            unsigned int tpc   = 0;
            unsigned int cstat = 0;
            theGeometry->PositionToTPC(pos, tpc, cstat);

            endT = nt;
            if (!foundStartPosition)
            {
                startT = endT;
                foundStartPosition = true;
            }   
        }
        catch(cet::exception &e){
            continue;
        }
    }
}

} // namespace lar_pandora
