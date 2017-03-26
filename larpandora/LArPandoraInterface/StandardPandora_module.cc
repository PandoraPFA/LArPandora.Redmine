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
   
    bool                m_uniqueInstanceSettings;   ///< Whether to enable unique configuration of each Pandora instance
    std::string         m_outputGeometryXmlFile;    ///< If provided, attempt to write collected geometry information to output xml file
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

#include <iostream>

namespace lar_pandora
{

StandardPandora::StandardPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
    m_uniqueInstanceSettings = pset.get<bool>("UniqueInstanceSettings", false);
    m_outputGeometryXmlFile = pset.get<std::string>("OutputGeometryXmlFile", "");
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
        throw cet::exception("LArPandora") << " Throwing exception - list of drift volumes is empty ";

    cet::search_path sp("FW_SEARCH_PATH");
    std::string fullConfigFileName;

    if (!sp.find_file(configFileName, fullConfigFileName))
        throw cet::exception("LArPandora") << " Failed to find xml configuration file " << configFileName << " in FW search path";

    const LArDriftVolume &driftVolume(m_driftVolumeList.front());

    m_pPrimaryPandora = this->CreateNewPandora();
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, fullConfigFileName));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora,
        new lar_content::LArPseudoLayerPlugin(driftVolume.GetWirePitchU(), driftVolume.GetWirePitchV(), driftVolume.GetWirePitchW())));

    // If only single drift volume, primary pandora instance will do all pattern recognition, rather than perform a particle stitching role
    if (1 == m_driftVolumeList.size())
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*m_pPrimaryPandora,
            new lar_content::LArRotationalTransformationPlugin(driftVolume.GetWireAngleU(), driftVolume.GetWireAngleV(), driftVolume.GetSigmaUVZ())));

        MultiPandoraApi::SetVolumeInfo(m_pPrimaryPandora, new VolumeInfo(0, "driftVolume",
            driftVolume.GetCenterX(), driftVolume.GetCenterY(), driftVolume.GetCenterZ(), 
	    driftVolume.GetWidthX(), driftVolume.GetWidthY(), driftVolume.GetWidthZ(), driftVolume.IsPositiveDrift()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreateDaughterPandoraInstances(const std::string &configFileName) const
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreateDaughterPandoraInstance(...) *** " << std::endl;

    if (!m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " Throwing exception - trying to create daughter Pandora instances in absence of primary instance ";

    for (const LArDriftVolume &driftVolume : m_driftVolumeList)
    {
        mf::LogDebug("LArPandora") << " Creating Pandora Daughter Instance: [" << driftVolume.GetVolumeID() << "]" << std::endl;

        std::ostringstream volumeIdString;
        volumeIdString << driftVolume.GetVolumeID();

        const pandora::Pandora *const pPandora = this->CreateNewPandora();
        MultiPandoraApi::AddDaughterPandoraInstance(m_pPrimaryPandora, pPandora);
        MultiPandoraApi::SetVolumeInfo(pPandora, new VolumeInfo(driftVolume.GetVolumeID(), "driftVolume_" + volumeIdString.str(),
            driftVolume.GetCenterX(), driftVolume.GetCenterY(), driftVolume.GetCenterZ(), 
            driftVolume.GetWidthX(), driftVolume.GetWidthY(), driftVolume.GetWidthZ(), driftVolume.IsPositiveDrift()));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora,
            new lar_content::LArPseudoLayerPlugin(driftVolume.GetWirePitchU(), driftVolume.GetWirePitchV(), driftVolume.GetWirePitchW())));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora,
            new lar_content::LArRotationalTransformationPlugin(driftVolume.GetWireAngleU(), driftVolume.GetWireAngleV(), driftVolume.GetSigmaUVZ())));

        std::string thisConfigFileName(configFileName);

        if (m_uniqueInstanceSettings)
        {
            const size_t insertPosition((thisConfigFileName.length() < 4) ? 0 : thisConfigFileName.length() - std::string(".xml").length());
            thisConfigFileName = thisConfigFileName.insert(insertPosition, volumeIdString.str());
        }

        cet::search_path sp("FW_SEARCH_PATH");
        std::string thisFullConfigFileName;

        if (!sp.find_file(thisConfigFileName, thisFullConfigFileName))
            throw cet::exception("LArPandora") << " Failed to find xml configuration file " << thisConfigFileName << " in FW search path";

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, thisFullConfigFileName));
    }
}

} // namespace lar_pandora
