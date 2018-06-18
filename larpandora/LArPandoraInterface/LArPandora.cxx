/**
 *  @file   larpandora/LArPandoraInterface/LArPandora.cxx
 *
 *  @brief  Base producer module for reconstructing recob::PFParticles from recob::Hits
 *
 */

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib/cpu_timer.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraInput.h"
#include "larpandora/LArPandoraInterface/LArPandoraOutput.h"

#include <iostream>
#include <limits>

namespace lar_pandora
{

LArPandora::LArPandora(fhicl::ParameterSet const &pset) :
    ILArPandora(pset),
    m_configFile(pset.get<std::string>("ConfigFile")),
    m_shouldRunAllHitsCosmicReco(pset.get<bool>("ShouldRunAllHitsCosmicReco")),
    m_shouldRunStitching(pset.get<bool>("ShouldRunStitching")),
    m_shouldRunCosmicHitRemoval(pset.get<bool>("ShouldRunCosmicHitRemoval")),
    m_shouldRunSlicing(pset.get<bool>("ShouldRunSlicing")),
    m_shouldRunNeutrinoRecoOption(pset.get<bool>("ShouldRunNeutrinoRecoOption")),
    m_shouldRunCosmicRecoOption(pset.get<bool>("ShouldRunCosmicRecoOption")),
    m_shouldPerformSliceId(pset.get<bool>("ShouldPerformSliceId")),
    m_printOverallRecoStatus(pset.get<bool>("PrintOverallRecoStatus", false)),
    m_generatorModuleLabel(pset.get<std::string>("GeneratorModuleLabel", "")),
    m_geantModuleLabel(pset.get<std::string>("GeantModuleLabel", "largeant")),
    m_hitfinderModuleLabel(pset.get<std::string>("HitFinderModuleLabel")),
    m_backtrackerModuleLabel(pset.get<std::string>("BackTrackerModuleLabel","")),
    m_enableProduction(pset.get<bool>("EnableProduction", true)),
    m_enableDetectorGaps(pset.get<bool>("EnableLineGaps", true)),
    m_enableMCParticles(pset.get<bool>("EnableMCParticles", false)),
    m_lineGapsCreated(false)
{
    m_inputSettings.m_useHitWidths = pset.get<bool>("UseHitWidths", true);
    m_inputSettings.m_useBirksCorrection = pset.get<bool>("UseBirksCorrection", false);
    m_inputSettings.m_uidOffset = pset.get<int>("UidOffset", 100000000);
    m_inputSettings.m_dx_cm = pset.get<double>("DefaultHitWidth", 0.5);
    m_inputSettings.m_int_cm = pset.get<double>("InteractionLength", 84.);
    m_inputSettings.m_rad_cm = pset.get<double>("RadiationLength", 14.);
    m_inputSettings.m_dEdX_mip = pset.get<double>("dEdXmip", 2.);
    m_inputSettings.m_mips_max = pset.get<double>("MipsMax", 50.);
    m_inputSettings.m_mips_if_negative = pset.get<double>("MipsIfNegative", 0.);
    m_inputSettings.m_mips_to_gev = pset.get<double>("MipsToGeV", 3.5e-4);
    m_inputSettings.m_recombination_factor = pset.get<double>("RecombinationFactor", 0.63);
    m_outputSettings.m_pProducer = this;
    m_outputSettings.m_shouldRunStitching = m_shouldRunStitching;

    if (m_enableProduction)
    {
        produces< std::vector<recob::PFParticle> >();
        produces< std::vector<recob::SpacePoint> >();
        produces< std::vector<recob::Cluster> >();
        produces< std::vector<recob::Vertex> >();
        produces< std::vector<larpandoraobj::PFParticleMetadata> >();

        produces< art::Assns<recob::PFParticle, larpandoraobj::PFParticleMetadata> >();
        produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
        produces< art::Assns<recob::PFParticle, recob::Cluster> >();
        produces< art::Assns<recob::PFParticle, recob::Vertex> >();
        produces< art::Assns<recob::SpacePoint, recob::Hit> >();
        produces< art::Assns<recob::Cluster, recob::Hit> >();

        if (m_outputSettings.m_shouldRunStitching)
        {
            produces< std::vector<anab::T0> >();
            produces< art::Assns<recob::PFParticle, anab::T0> >();
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::beginJob()
{
    LArDriftVolumeList driftVolumeList;
    LArPandoraGeometry::LoadGeometry(driftVolumeList, m_driftVolumeMap);

    this->CreatePandoraInstances();

    if (!m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " LArPandora::beginJob - failed to create primary Pandora instance " << std::endl;

    m_inputSettings.m_pPrimaryPandora = m_pPrimaryPandora;
    m_outputSettings.m_pPrimaryPandora = m_pPrimaryPandora;

    // Pass basic LArTPC information to pandora instances
    LArPandoraInput::CreatePandoraLArTPCs(m_inputSettings, driftVolumeList);

    // If using global drift volume approach, pass details of gaps between daughter volumes to the pandora instance
    if (m_enableDetectorGaps)
    {
        LArDetectorGapList listOfGaps;
        LArPandoraGeometry::LoadDetectorGaps(listOfGaps);
        LArPandoraInput::CreatePandoraDetectorGaps(m_inputSettings, driftVolumeList, listOfGaps);
    }

    // Parse Pandora settings xml files
    this->ConfigurePandoraInstances();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::produce(art::Event &evt)
{
    IdToHitMap idToHitMap;
    this->CreatePandoraInput(evt, idToHitMap);
    this->RunPandoraInstances();
    this->ProcessPandoraOutput(evt, idToHitMap);
    this->ResetPandoraInstances();
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
    RawMCParticleVector generatorArtMCParticleVector;
    MCTruthToMCParticles artMCTruthToMCParticles;
    MCParticlesToMCTruth artMCParticlesToMCTruth;

    LArPandoraHelper::CollectHits(evt, m_hitfinderModuleLabel, artHits);

    if (m_enableMCParticles && !evt.isRealData())
    {
        LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCParticleVector);

        if (!m_generatorModuleLabel.empty())
            LArPandoraHelper::CollectGeneratorMCParticles(evt, m_generatorModuleLabel, generatorArtMCParticleVector);

        LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth);

        LArPandoraHelper::CollectSimChannels(evt, m_geantModuleLabel, artSimChannels);
        if (!artSimChannels.empty())
        {
            LArPandoraHelper::BuildMCParticleHitMaps(artHits, artSimChannels, artHitsToTrackIDEs);
        }
        else
        {
            if (m_backtrackerModuleLabel.empty())
              throw cet::exception("LArPandora") << " LArPandora::CreatePandoraInput - no sim channels found, backtracker module must be set in FHiCL " << std::endl;

            LArPandoraHelper::BuildMCParticleHitMaps(evt, m_hitfinderModuleLabel, m_backtrackerModuleLabel, artHitsToTrackIDEs);
        }
    }

    LArPandoraInput::CreatePandoraHits2D(m_inputSettings, m_driftVolumeMap, artHits, idToHitMap);

    if (m_enableMCParticles && !evt.isRealData())
    {
        LArPandoraInput::CreatePandoraMCParticles(m_inputSettings, artMCTruthToMCParticles, artMCParticlesToMCTruth, generatorArtMCParticleVector);
        LArPandoraInput::CreatePandoraMCLinks2D(m_inputSettings, idToHitMap, artHitsToTrackIDEs);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandora::ProcessPandoraOutput(art::Event &evt, const IdToHitMap &idToHitMap)
{
    if (m_enableProduction)
        LArPandoraOutput::ProduceArtOutput(m_outputSettings, idToHitMap, evt);
}

} // namespace lar_pandora
