/**
 *  @file   larpandora/LArPandoraInterface/LArPandora.cxx
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/cpu_timer.h"

#include "lardata/ArtDataHelper/MVAReader.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include <iostream>

namespace lar_pandora
{

LArPandora::LArPandora(fhicl::ParameterSet const &pset) :
    ILArPandora(pset)
{
    m_configFile = pset.get<std::string>("ConfigFile");
    m_stitchingConfigFile = pset.get<std::string>("StitchingConfigFile", "No_File_Provided");

    // Configuration for this module
    m_runStitchingInstance = pset.get<bool>("RunStitchingInstance", false);
    m_enableProduction = pset.get<bool>("EnableProduction", true);
    m_enableDetectorGaps = pset.get<bool>("EnableLineGaps", true);
    m_enableMCParticles = pset.get<bool>("EnableMCParticles", false);

    // Input labels for this module
    m_geantModuleLabel = pset.get<std::string>("GeantModuleLabel", "largeant");
    m_hitfinderModuleLabel = pset.get<std::string>("HitFinderModuleLabel");
    m_mvaModuleLabel = pset.get<std::string>("MVAModuleLabel", "");

    // Settings for LArPandoraGeometry
    m_geometrySettings.m_globalCoordinates = pset.get<bool>("UseGlobalCoordinates", false);  // Keep separate drift volumes but interchange U and V
    m_geometrySettings.m_globalDriftVolume = pset.get<bool>("UseGlobalDriftVolume", false);  // Transform to a single global drift volume

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
    m_outputSettings.m_buildStitchedParticles = (m_runStitchingInstance && pset.get<bool>("BuildStitchedParticles", false));

    // Status flags for this module
    m_lineGapsCreated = false;

    if (m_enableProduction)
    {
        produces< std::vector<recob::PFParticle> >();
        produces< std::vector<recob::SpacePoint> >();
        produces< std::vector<recob::Cluster> >();
        produces< std::vector<recob::Vertex> >();

        produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
        produces< art::Assns<recob::PFParticle, recob::Cluster> >();
        produces< art::Assns<recob::PFParticle, recob::Vertex> >();
        produces< art::Assns<recob::SpacePoint, recob::Hit> >();
        produces< art::Assns<recob::Cluster, recob::Hit> >();

        if (m_outputSettings.m_buildStitchedParticles)
        {
            produces< std::vector<anab::T0> >();
            produces< art::Assns<recob::PFParticle, anab::T0> >();
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
    // Load geometry to determine number of required Pandora instances
    LArPandoraGeometry::LoadGeometry(m_geometrySettings, m_driftVolumeList, m_driftVolumeMap);

    this->CreatePandoraInstances();

    if (!m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " LArPandora::beginJob - failed to create primary Pandora instance " << std::endl;

    m_inputSettings.m_pPrimaryPandora = m_pPrimaryPandora;
    m_outputSettings.m_pPrimaryPandora = m_pPrimaryPandora;

    // Pass basic LArTPC information to pandora instances
    LArPandoraInput::CreatePandoraLArTPCs(m_inputSettings, m_driftVolumeList);

    // If using global drift volume approach, pass details of gaps between daughter volumes to the pandora instance
    if (m_enableDetectorGaps && m_geometrySettings.m_globalDriftVolume)
    {
        LArDetectorGapList listOfGaps;
        LArPandoraGeometry::LoadDetectorGaps(m_geometrySettings, listOfGaps);
        LArPandoraInput::CreatePandoraDetectorGaps(m_inputSettings, m_driftVolumeList, listOfGaps);
    }

    // Parse Pandora settings xml files
    this->ConfigurePandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::produce(art::Event &evt)
{
    mf::LogInfo("LArPandora") << " *** LArPandora::produce(...)  [Run=" << evt.run() << ", Event=" << evt.id().event() << "] *** " << std::endl;

    IdToHitMap idToHitMap;
    this->CreatePandoraInput(evt, idToHitMap);
    this->RunPandoraInstances();
    this->ProcessPandoraOutput(evt, idToHitMap);
    this->ResetPandoraInstances();

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

    std::unique_ptr<anab::MVAReader<recob::Hit, 4> > pHitResults;

    if (!m_mvaModuleLabel.empty())
    {
        pHitResults = anab::MVAReader<recob::Hit, 4>::create(evt, m_mvaModuleLabel);

        if (!pHitResults)
            mf::LogWarning("LArPandora") << "No MVA results available." << std::endl;

        if (pHitResults && ((pHitResults->getIndex("em") < 0) || (pHitResults->getIndex("track") < 0) || (pHitResults->getIndex("michel") < 0) || (pHitResults->getIndex("none") < 0)))
            throw cet::exception("LArPandora") << "MVA Data present, but no em/track/michel/none-labelled columns in MVA data products." << std::endl;
    }

    LArPandoraInput::CreatePandoraHits2D(m_inputSettings, m_driftVolumeMap, artHits, idToHitMap, pHitResults);

    if (m_enableMCParticles && !evt.isRealData())
    {
        LArPandoraInput::CreatePandoraMCParticles(m_inputSettings, m_driftVolumeMap, artMCTruthToMCParticles, artMCParticlesToMCTruth);
        LArPandoraInput::CreatePandoraMCLinks2D(m_inputSettings, m_driftVolumeMap, idToHitMap, artHitsToTrackIDEs);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ProcessPandoraOutput(art::Event &evt, const IdToHitMap &idToHitMap)
{
    if (m_enableProduction)
        LArPandoraOutput::ProduceArtOutput(m_outputSettings, idToHitMap, evt);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::RunPandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** LArPandora::RunPandoraInstances() *** " << std::endl;

    const PandoraInstanceList &daughterInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(m_pPrimaryPandora));

    for (const pandora::Pandora *const pPandora : daughterInstances)
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPandora));

    if (m_runStitchingInstance || daughterInstances.empty())
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pPrimaryPandora));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ResetPandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** LArPandora::ResetPandoraInstances() *** " << std::endl;

    const PandoraInstanceList &daughterInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(m_pPrimaryPandora));

    for (const pandora::Pandora *const pPandora : daughterInstances)
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPandora));

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

    return pPandora;
}

} // namespace lar_pandora
