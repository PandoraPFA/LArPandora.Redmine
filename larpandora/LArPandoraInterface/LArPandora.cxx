/**
 *  @file   larpandora/LArPandoraInterface/LArPandora.cxx
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "TTree.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include "larpandora/LArPandoraShowers/PCAShowerParticleBuildingAlgorithm.h"

#include <iostream>

namespace lar_pandora
{

LArPandora::LArPandora(fhicl::ParameterSet const &pset) :
    ILArPandora(pset)
{
    m_configFile = pset.get<std::string>("ConfigFile");
    m_stitchingConfigFile = pset.get<std::string>("StitchingConfigFile", "No_File_Provided");

    // prepare the optional cluster energy algorithm
    if (pset.has_key("ShowerEnergy") && pset.is_key_to_table("ShowerEnergy"))
    {
        m_showerEnergyAlg = std::make_unique<calo::LinearEnergyAlg>(pset.get<fhicl::ParameterSet>("ShowerEnergy"));
    }
    else mf::LogWarning("LArPandora") << "No shower energy calibration set up.";

    // Configuration for this module
    m_runStitchingInstance = pset.get<bool>("RunStitchingInstance", false);
    m_enableProduction = pset.get<bool>("EnableProduction", true);
    m_enableDetectorGaps = pset.get<bool>("EnableLineGaps", true);
    m_enableMCParticles = pset.get<bool>("EnableMCParticles", false);
    m_enableMonitoring = pset.get<bool>("EnableMonitoring", false);

    // Input labels for this module
    m_geantModuleLabel = pset.get<std::string>("GeantModuleLabel", "largeant");
    m_hitfinderModuleLabel = pset.get<std::string>("HitFinderModuleLabel", "gaushit");
    m_pandoraModuleLabel = pset.get<std::string>("PFParticleModuleLabel", "pandora");

    // Settings for LArPandoraGeometry
    m_geometrySettings.m_globalCoordinates = pset.get<bool>("UseGlobalCoordinates", false);  // Keep separate drift volumes but interchange U and V
    m_geometrySettings.m_globalDriftVolume = pset.get<bool>("UseGlobalDriftVolume", false);  // Transform to a single global drift volume
    m_geometrySettings.m_printGeometry = pset.get<bool>("PrintGeometry", false);

    // Settings for LArPandoraInput
    m_inputSettings.m_useHitWidths = pset.get<bool>("UseHitWidths", true);
    m_inputSettings.m_uidOffset = pset.get<int>("UidOffset", 100000000);
    m_inputSettings.m_dx_cm = pset.get<double>("DefaultHitWidth", 0.5);
    m_inputSettings.m_int_cm = pset.get<double>("InteractionLength", 84.0);
    m_inputSettings.m_rad_cm = pset.get<double>("RadiationLength", 14.0);
    m_inputSettings.m_dEdX_max = pset.get<double>("dEdXmax", 25.0);
    m_inputSettings.m_dEdX_mip = pset.get<double>("dEdXmip", 2.0);
    m_inputSettings.m_mips_to_gev = pset.get<double>("MipsToGeV", 3.5e-4);
    m_inputSettings.m_recombination_factor = pset.get<double>("RecombinationFactor", 0.63);
    m_inputSettings.m_globalViews = (m_geometrySettings.m_globalCoordinates || m_geometrySettings.m_globalDriftVolume);
    m_inputSettings.m_truncateReadout = m_geometrySettings.m_globalDriftVolume;

    // Settings for LArPandoraOutput
    m_outputSettings.m_pProducer = this;
    m_outputSettings.m_buildTracks = pset.get<bool>("BuildTracks", true);
    m_outputSettings.m_buildShowers = pset.get<bool>("BuildShowers", true);
    m_outputSettings.m_buildStitchedParticles = (m_runStitchingInstance && pset.get<bool>("BuildStitchedParticles", false));
    m_outputSettings.m_showerEnergyAlg = m_showerEnergyAlg.get(); // may be nullptr

    // Status flags for this module
    m_lineGapsCreated = false;

    if (m_enableProduction)
    {
        produces< std::vector<recob::PFParticle> >();
        produces< std::vector<recob::SpacePoint> >();
        produces< std::vector<recob::Cluster> >();
        produces< std::vector<recob::Seed> >();
        produces< std::vector<recob::Vertex> >();

        produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
        produces< art::Assns<recob::PFParticle, recob::Cluster> >();
        produces< art::Assns<recob::PFParticle, recob::Seed> >();
        produces< art::Assns<recob::PFParticle, recob::Vertex> >();
        produces< art::Assns<recob::SpacePoint, recob::Hit> >();
        produces< art::Assns<recob::Cluster, recob::Hit> >();
        produces< art::Assns<recob::Seed, recob::Hit> >();

        if (m_outputSettings.m_buildTracks)
        {
            produces< std::vector<recob::Track> >();
            produces< art::Assns<recob::PFParticle, recob::Track> >();
            produces< art::Assns<recob::Track, recob::Hit> >();

            if (m_outputSettings.m_buildStitchedParticles)
            {
                produces< std::vector<anab::T0> >();
                produces< art::Assns<recob::Track, anab::T0> >();
            }
        }

        if (m_outputSettings.m_buildShowers)
        {
            produces< std::vector<recob::Shower> >();
            produces< std::vector<recob::PCAxis> >();
            produces< art::Assns<recob::PFParticle, recob::Shower> >();
            produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
            produces< art::Assns<recob::Shower, recob::Hit> >();
            produces< art::Assns<recob::Shower, recob::PCAxis> >();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPandora::~LArPandora()
{
    this->DeletePandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::beginJob()
{
    if (m_enableMonitoring)
        this->InitializeMonitoring();

    // Load geometry and then create Pandora instances
    LArPandoraGeometry::LoadGeometry(m_geometrySettings, m_driftVolumeList, m_driftVolumeMap);

    this->CreatePandoraInstances();

    if (!m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " LArPandora::beginJob - failed to create primary Pandora instance " << std::endl;

    m_inputSettings.m_pPrimaryPandora = m_pPrimaryPandora;
    m_outputSettings.m_pPrimaryPandora = m_pPrimaryPandora;

    // Load gaps associated with dead regions between drift volumes
    if (m_enableDetectorGaps)
    {
        LArDetectorGapList listOfGaps;
        LArPandoraGeometry::LoadDetectorGaps(m_geometrySettings, listOfGaps);
        LArPandoraInput::CreatePandoraDetectorGaps(m_inputSettings, listOfGaps);
    }

    // Print the configuration of the algorithm at the beginning of the job;
    // the algorithm does not need to be set up for this.
    if (m_showerEnergyAlg) {
      mf::LogInfo log("LArPandora");
      log << "Energy shower settings: ";
      m_showerEnergyAlg->DumpConfiguration(log, "  ", "");
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::produce(art::Event &evt)
{
    mf::LogInfo("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] *** " << std::endl;

    // we set up the algorithm on each new event, in case the services have changed:
    if (m_showerEnergyAlg) {
      m_showerEnergyAlg->setup(
        *(lar::providerFrom<detinfo::DetectorPropertiesService>()),
        *(lar::providerFrom<detinfo::DetectorClocksService>()),
        *(lar::providerFrom<geo::Geometry>())
        );
    } // if


    IdToHitMap idToHitMap;
    this->CreatePandoraInput(evt, idToHitMap);
    this->RunPandoraInstances();
    this->ProcessPandoraOutput(evt, idToHitMap);
    this->ResetPandoraInstances();

    if (m_enableMonitoring)
    {
        m_run = evt.run();
        m_event = evt.id().event();
        m_pRecoTree->Fill();
    }

    mf::LogInfo("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] DONE! *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::CreatePandoraInput(art::Event &evt, IdToHitMap &idToHitMap)
{
    // ATTN Should complete gap creation in begin job callback, but channel status service functionality unavailable at that point
    if (!m_lineGapsCreated && m_enableDetectorGaps)
    {
        LArPandoraInput::CreatePandoraReadoutGaps(m_inputSettings, m_driftVolumeMap);
        m_lineGapsCreated = true;
    }

    cet::cpu_timer theClock;

    if (m_enableMonitoring)
        theClock.start();

    HitVector artHits;
    SimChannelVector artSimChannels;
    HitsToTrackIDEs artHitsToTrackIDEs;
    MCParticleVector artMCParticleVector;
    MCTruthToMCParticles artMCTruthToMCParticles;
    MCParticlesToMCTruth artMCParticlesToMCTruth;

    LArPandoraHelper::CollectHits(evt, m_hitfinderModuleLabel, artHits);

    if (m_enableMCParticles && !evt.isRealData())
    {
        LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCParticleVector);
        LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth);
        LArPandoraHelper::CollectSimChannels(evt, m_geantModuleLabel, artSimChannels);
        LArPandoraHelper::BuildMCParticleHitMaps(artHits, artSimChannels, artHitsToTrackIDEs);
    }

    if (m_enableMonitoring)
    {
        theClock.stop();
        m_collectionTime = theClock.accumulated_real_time();
        theClock.reset();
        theClock.start();
    }

    LArPandoraInput::CreatePandoraHits2D(m_inputSettings, m_driftVolumeMap, artHits, idToHitMap);

    if (m_enableMCParticles && !evt.isRealData())
    {
        LArPandoraInput::CreatePandoraMCParticles(m_inputSettings, m_driftVolumeMap, artMCTruthToMCParticles, artMCParticlesToMCTruth);
        LArPandoraInput::CreatePandoraMCParticles2D(m_inputSettings, m_driftVolumeMap, artMCParticleVector);
        LArPandoraInput::CreatePandoraMCLinks2D(m_inputSettings, m_driftVolumeMap, idToHitMap, artHitsToTrackIDEs);
    }

    if (m_enableMonitoring)
    {
        theClock.stop();
        m_inputTime = theClock.accumulated_real_time();
        m_hits = static_cast<int>(artHits.size());
        m_pandoraHits = static_cast<int>(idToHitMap.size());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ProcessPandoraOutput(art::Event &evt, const IdToHitMap &idToHitMap)
{
    cet::cpu_timer theClock;

    if (m_enableMonitoring)
        theClock.start();

    if (m_enableProduction)
        LArPandoraOutput::ProduceArtOutput(m_outputSettings, idToHitMap, evt);

    if (m_enableMonitoring)
    {
        theClock.stop();
        m_outputTime = theClock.accumulated_real_time();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::RunPandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** LArPandora::RunPandoraInstances() *** " << std::endl;
    cet::cpu_timer theClock;

    if (m_enableMonitoring)
        theClock.start();

    const PandoraInstanceList &daughterInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(m_pPrimaryPandora));

    for (const pandora::Pandora *const pPandora : daughterInstances)
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPandora));
        this->SetParticleX0Values(pPandora);
    }

    if (m_runStitchingInstance || daughterInstances.empty())
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPrimaryPandora));

    if (m_enableMonitoring)
    {
        theClock.stop();
        m_processTime = theClock.accumulated_real_time();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ResetPandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ResetPandoraInstances() *** " << std::endl;

    const PandoraInstanceList &daughterInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(m_pPrimaryPandora));

    for (const pandora::Pandora *const pPandora : daughterInstances)
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPandora));
        MultiPandoraApi::ClearParticleX0Map(pPandora);
    }

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pPrimaryPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::DeletePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** LArPandora::DeletePandoraInstances() *** " << std::endl;

    if (m_pPrimaryPandora)
      MultiPandoraApi::DeletePandoraInstances(m_pPrimaryPandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Pandora *LArPandora::CreateNewPandora() const
{
    const pandora::Pandora *const pPandora = new pandora::Pandora();
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*pPandora,
        "LArPCAShowerParticleBuilding", new lar_pandora_showers::PCAShowerParticleBuildingAlgorithm::Factory));

    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::SetParticleX0Values(const pandora::Pandora *const pPandora) const
{
    // ATTN Just a placeholder for a proper treatment
    const pandora::PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pPandora, pPfoList));

    pandora::PfoList connectedPfoList;
    lar_content::LArPfoHelper::GetAllConnectedPfos(*pPfoList, connectedPfoList);

    for (const pandora::ParticleFlowObject *const pPfo : connectedPfoList)
        MultiPandoraApi::SetParticleX0(pPandora, pPfo, 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::InitializeMonitoring()
{
    art::ServiceHandle<art::TFileService> tfs;
    m_pRecoTree = tfs->make<TTree>("monitoring", "LAr Reco");
    m_pRecoTree->Branch("run", &m_run, "run/I");
    m_pRecoTree->Branch("event", &m_event, "event/I");
    m_pRecoTree->Branch("hits", &m_hits, "hits/I");
    m_pRecoTree->Branch("pandoraHits", &m_pandoraHits, "pandoraHits/I");
    m_pRecoTree->Branch("collectionTime", &m_collectionTime, "collectionTime/F");
    m_pRecoTree->Branch("inputTime", &m_inputTime, "inputTime/F");
    m_pRecoTree->Branch("processTime", &m_processTime, "processTime/F");
    m_pRecoTree->Branch("outputTime", &m_outputTime, "outputTime/F");
}

} // namespace lar_pandora
