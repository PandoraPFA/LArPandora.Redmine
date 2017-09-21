/**
 *  @file   larpandora/LArPandoraInterface/StandardPandora_module.cc
 *
 *  @brief  A generic LArPandora ART Producer module intended to work on ALL LAr-TPC wire-readout experiments.
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"

#include <string>

namespace lar_pandora
{

/**
 *  @brief  StandardPandora class
 */
class StandardPandora : public LArPandora
{
public:

    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    StandardPandora(fhicl::ParameterSet const &pset);

private:
    void CreatePandoraInstances();

    /**
     *  @brief  Create primary pandora instance
     *
     *  @param  configFileName the pandora settings config file name
     */
    void CreatePrimaryPandoraInstance(const std::string &configFileName);

    /**
     *  @brief  Create daughter pandora instances
     *
     *  @param  configFileName the pandora settings config file name
     */
    void CreateDaughterPandoraInstances(const std::string &configFileName) const;

    /**
     *  @brief  Pass external steering parameters, read from fhicl parameter set, to LArParent Pandora algorithm
     * 
     *  @param  pPandora the address of the relevant pandora instance
     */
    void ProvideExternalSteeringParameters(const pandora::Pandora *const pPandora) const;

    bool                m_uniqueInstanceSettings;           ///< Whether to enable unique configuration of each Pandora instance
    std::string         m_outputGeometryXmlFile;            ///< If provided, attempt to write collected geometry information to output xml file

    bool                m_shouldRunAllHitsCosmicReco;       ///< Steering: whether to run all hits cosmic-ray reconstruction
    bool                m_shouldRunCosmicHitRemoval;        ///< Steering: whether to remove hits from tagged cosmic-rays
    bool                m_shouldRunSlicing;                 ///< Steering: whether to slice events into separate regions for processing
    bool                m_shouldRunNeutrinoRecoOption;      ///< Steering: whether to run neutrino reconstruction for each slice
    bool                m_shouldRunCosmicRecoOption;        ///< Steering: whether to run cosmic-ray reconstruction for each slice
    bool                m_shouldIdentifyNeutrinoSlice;      ///< Steering: whether to identify most appropriate neutrino slice
    bool                m_printOverallRecoStatus;           ///< Steering: whether to print current operation status messages
};

DEFINE_ART_MODULE(StandardPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib/exception.h"

#include "larcore/Geometry/Geometry.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"
#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArUtility/ParentAlgorithm.h"

#include <iostream>

namespace lar_pandora
{

StandardPandora::StandardPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
    m_uniqueInstanceSettings = pset.get<bool>("UniqueInstanceSettings", false);
    m_outputGeometryXmlFile = pset.get<std::string>("OutputGeometryXmlFile", "");

    m_shouldRunAllHitsCosmicReco = pset.get<bool>("ShouldRunAllHitsCosmicReco");
    m_shouldRunCosmicHitRemoval = pset.get<bool>("ShouldRunCosmicHitRemoval");
    m_shouldRunSlicing = pset.get<bool>("ShouldRunSlicing");
    m_shouldRunNeutrinoRecoOption = pset.get<bool>("ShouldRunNeutrinoRecoOption");
    m_shouldRunCosmicRecoOption = pset.get<bool>("ShouldRunCosmicRecoOption");
    m_shouldIdentifyNeutrinoSlice = pset.get<bool>("ShouldIdentifyNeutrinoSlice");
    m_printOverallRecoStatus = pset.get<bool>("PrintOverallRecoStatus", false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreatePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreatePandoraInstances(...) [BEGIN] *** " << std::endl;

    // For multiple drift volumes, create a Pandora instance for each drift volume and an additional instance for stitching drift volumes.
    // For single drift volumes, just create a single Pandora instance
    if (m_driftVolumeList.size() > 1)
    {
        this->CreatePrimaryPandoraInstance(m_stitchingConfigFile);
        this->CreateDaughterPandoraInstances(m_configFile);
    }
    else
    {
        this->CreatePrimaryPandoraInstance(m_configFile);
    }

    // Write out geometry here
    if (!m_outputGeometryXmlFile.empty())
        LArPandoraGeometry::WriteGeometry(m_outputGeometryXmlFile, m_driftVolumeList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreatePrimaryPandoraInstance(const std::string &configFileName)
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreatePrimaryPandoraInstance(...) *** " << std::endl;

    if (m_driftVolumeList.empty())
        throw cet::exception("LArPandora") << " StandardPandora::CreatePrimaryPandoraInstance --- list of drift volumes is empty ";

    m_pPrimaryPandora = this->CreateNewPandora();
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);

    // If only single drift volume, primary pandora instance will do all pattern recognition, rather than perform a particle stitching role
    if (1 == m_driftVolumeList.size())
    {
        this->ProvideExternalSteeringParameters(m_pPrimaryPandora);
        MultiPandoraApi::SetVolumeId(m_pPrimaryPandora, 0);
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*m_pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*m_pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));
    }

    cet::search_path sp("FW_SEARCH_PATH");
    std::string fullConfigFileName;

    if (!sp.find_file(configFileName, fullConfigFileName))
        throw cet::exception("LArPandora") << " StandardPandora::CreatePrimaryPandoraInstance --- Failed to find xml configuration file " << configFileName << " in FW search path";

    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, fullConfigFileName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreateDaughterPandoraInstances(const std::string &configFileName) const
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreateDaughterPandoraInstance(...) *** " << std::endl;

    if (!m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " StandardPandora::CreateDaughterPandoraInstances --- trying to create daughter Pandora instances in absence of primary instance ";

    if (m_driftVolumeList.size() < 2)
        return;

    for (const LArDriftVolume &driftVolume : m_driftVolumeList)
    {
        const pandora::Pandora *const pPandora = this->CreateNewPandora();
        MultiPandoraApi::AddDaughterPandoraInstance(m_pPrimaryPandora, pPandora);
        MultiPandoraApi::SetVolumeId(pPandora, driftVolume.GetVolumeID());
        this->ProvideExternalSteeringParameters(pPandora);

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));

        std::string thisConfigFileName(configFileName);

        if (m_uniqueInstanceSettings)
        {
            std::ostringstream volumeIdString;
            volumeIdString << driftVolume.GetVolumeID();
            const size_t insertPosition((thisConfigFileName.length() < 4) ? 0 : thisConfigFileName.length() - std::string(".xml").length());
            thisConfigFileName = thisConfigFileName.insert(insertPosition, volumeIdString.str());
        }

        cet::search_path sp("FW_SEARCH_PATH");
        std::string thisFullConfigFileName;

        if (!sp.find_file(thisConfigFileName, thisFullConfigFileName))
            throw cet::exception("LArPandora") << "  StandardPandora::CreateDaughterPandoraInstances --- failed to find xml configuration file " << thisConfigFileName << " in FW search path";

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, thisFullConfigFileName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::ProvideExternalSteeringParameters(const pandora::Pandora *const pPandora) const
{
    auto *const pEventSteeringParameters = new lar_content::ParentAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldIdentifyNeutrinoSlice = m_shouldIdentifyNeutrinoSlice;
    pEventSteeringParameters->m_printOverallRecoStatus = m_printOverallRecoStatus;
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::ExternallyConfiguredAlgorithm::SetExternalParameters(*pPandora, "LArParent", pEventSteeringParameters));
}

} // namespace lar_pandora
