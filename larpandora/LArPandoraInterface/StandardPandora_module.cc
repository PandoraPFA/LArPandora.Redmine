/**
 *  @file   larpandora/LArPandoraInterface/StandardPandora_module.cc
 *
 *  @brief  A generic LArPandora ART Producer module intended to work on ALL LAr-TPC wire-readout experiments.
 *          This single module is intended to take the place of the current suite of experiment-specific modules.
 *
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

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
    int GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const;

    /**
     *  @brief Output a unique ID number for each TPC
     *
     *  @param cryostat ID number
     *  @param tpc ID number
     */
    unsigned int GetTpcIdNumber(const unsigned int cryostat, const unsigned int tpc) const;

    /**
     *  @brief Load the geometry information needed to run the Pandora reconstruction
     */
    void LoadGeometry();

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
    void CreateDaughterPandoraInstances(const std::string &configFileName);

    LArDriftVolumeList  m_driftVolumeList;   // list of drift volumes for this geometry
    LArDriftVolumeMap   m_driftVolumeMap;    // mapping from tpcVolumeID to driftVolumeID

    bool     m_printGeometry;
    bool     m_useShortVolume;   // Historical DUNE 35t config parameter - use short drift volume (positive drift)
    bool     m_useLongVolume;    // Historical DUNE 35t config parameter - use long drift volume (negative drift)
    bool     m_useLeftVolume;    // Historical DUNE 4-APA config parameter - use left drift volume (negative drift)
    bool     m_useRightVolume;   // Historical DUNE 4-APA config parameter - use right drift volume (positive drift)
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

namespace lar_pandora
{

StandardPandora::StandardPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
    m_printGeometry = pset.get<bool>("PrintGeometry", false);
    m_useShortVolume = pset.get<bool>("UseShortVolume", true);
    m_useLongVolume = pset.get<bool>("UseLongVolume", true);
    m_useLeftVolume = pset.get<bool>("UseLeftVolume", true);
    m_useRightVolume = pset.get<bool>("UseRightVolume", true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int StandardPandora::GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const
{
    LArDriftVolumeMap::const_iterator iter = m_driftVolumeMap.find(this->GetTpcIdNumber(cryostat, tpc));

    if (m_driftVolumeMap.end() == iter)
        throw cet::exception("LArPandora") << " Throwing exception - found a TPC that doesn't belong to a drift volume";
 
    const LArDriftVolume &driftVolume = iter->second;
    return driftVolume.GetVolumeID();
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int StandardPandora::GetTpcIdNumber(const unsigned int cryostat, const unsigned int tpc) const
{
    return 10000 * cryostat + tpc;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::LoadGeometry()
{
    if (!m_driftVolumeList.empty() || !m_driftVolumeMap.empty())
        throw cet::exception("LArPandora") << " Throwing exception - list of drift volumes already exists ";

    LArPandoraGeometry::LoadGeometry(m_driftVolumeList);
 
    if (m_printGeometry)
        LArPandoraGeometry::PrintGeometry(m_driftVolumeList);

    for (LArDriftVolumeList::const_iterator iter1 = m_driftVolumeList.begin(), iterEnd1 = m_driftVolumeList.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArDriftVolume &driftVolume = *iter1;
        const LArTpcVolumeList &tpcVolumeList = driftVolume.GetTpcVolumeList();

        for (LArTpcVolumeList::const_iterator iter2 = tpcVolumeList.begin(), iterEnd2 = tpcVolumeList.end(); iter2 != iterEnd2; ++iter2)
	{
	    const LArTpcVolume &tpcVolume = *iter2;
	    m_driftVolumeMap.insert(LArDriftVolumeMap::value_type(this->GetTpcIdNumber(tpcVolume.GetCryostat(), tpcVolume.GetTpc()), driftVolume));
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreatePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreatePandoraInstances(...) [BEGIN] *** " << std::endl;
    
    // Load (and print) geometry
    this->LoadGeometry();

    // Ensure that xml configuration files are available
    const bool isMultiDrift(m_driftVolumeList.size() > 1);

    cet::search_path sp("FW_SEARCH_PATH");
    std::string stitchingConfigFileName, configFileName;

    if (!sp.find_file(m_configFile, configFileName))
    {
        throw cet::exception("LArPandora") << " Failed to find xml configuration file " << m_configFile << " in FW seach path";
    }

    if (isMultiDrift && !sp.find_file(m_stitchingConfigFile, stitchingConfigFileName))
    {
        throw cet::exception("LArPandora") << " Failed to find xml configuration file " << m_stitchingConfigFile << " in FW search path";
    }

    // For multiple drift volumes, create a Pandora instance for each drift volume and an additional instance for stitching drift volumes.
    if (isMultiDrift)
    {
        this->CreatePrimaryPandoraInstance(stitchingConfigFileName);
        this->CreateDaughterPandoraInstances(configFileName);
    }

    // For single drift volumes, just create a single Pandora instance
    else
    {
        this->CreatePrimaryPandoraInstance(configFileName);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreatePrimaryPandoraInstance(const std::string &configFileName)
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreatePrimaryPandoraInstance(...) *** " << std::endl;

    if (m_driftVolumeList.empty())
        throw cet::exception("LArPandora") << " Throwing exception - list of drift volumes is empty ";

    const LArDriftVolume &driftVolume = *(m_driftVolumeList.begin());

    // Common settings if detector has single and multiple drift volumes
    m_pPrimaryPandora = this->CreateNewPandora();
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, configFileName));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora,
        new lar_content::LArPseudoLayerPlugin(driftVolume.GetWirePitchU(), driftVolume.GetWirePitchV(), driftVolume.GetWirePitchW())));

    // Additional settings if detector has just a single drift volume
    if (1 == m_driftVolumeList.size())
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*m_pPrimaryPandora,
	    new lar_content::LArRotationalTransformationPlugin(driftVolume.GetWireAngleU(), driftVolume.GetWireAngleV(), driftVolume.GetSigmaUVZ())));
        MultiPandoraApi::SetVolumeInfo(m_pPrimaryPandora, new VolumeInfo(0, "driftVolume",
	    pandora::CartesianVector(driftVolume.GetCenterX(), driftVolume.GetCenterY(), driftVolume.GetCenterZ()), driftVolume.IsPositiveDrift()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void StandardPandora::CreateDaughterPandoraInstances(const std::string &configFileName)
{
    mf::LogDebug("LArPandora") << " *** StandardPandora::CreateDaughterPandoraInstance(...) *** " << std::endl;

    if (!m_pPrimaryPandora)
        throw cet::exception("LArPandora") << " Throwing exception - trying to create daughter Pandora instances in absence of primary instance ";

    for (LArDriftVolumeList::const_iterator iter = m_driftVolumeList.begin(), iterEnd = m_driftVolumeList.end(); iter != iterEnd; ++iter)
    {
        const LArDriftVolume &driftVolume = *iter;

        // Check historical DUNE config parameters
        if ((2 == m_driftVolumeList.size()) &&
            ((true == driftVolume.IsPositiveDrift() && (false == m_useShortVolume || false == m_useRightVolume)) ||
             (false == driftVolume.IsPositiveDrift() && (false == m_useLongVolume || false == m_useLeftVolume))))
	    continue;

	mf::LogDebug("LArPandora") << " Creating Pandora Daughter Instance: [" << driftVolume.GetVolumeID() << "]" << std::endl;

        std::ostringstream volumeIdString("driftVolume_");
        volumeIdString << driftVolume.GetVolumeID();

        const pandora::Pandora *const pPandora = this->CreateNewPandora();
        MultiPandoraApi::AddDaughterPandoraInstance(m_pPrimaryPandora, pPandora);
        MultiPandoraApi::SetVolumeInfo(pPandora, new VolumeInfo(driftVolume.GetVolumeID(), volumeIdString.str(),
	    pandora::CartesianVector(driftVolume.GetCenterX(), driftVolume.GetCenterY(), driftVolume.GetCenterZ()), driftVolume.IsPositiveDrift()));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora,
	    new lar_content::LArPseudoLayerPlugin(driftVolume.GetWirePitchU(), driftVolume.GetWirePitchV(), driftVolume.GetWirePitchW())));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora,
	    new lar_content::LArRotationalTransformationPlugin(driftVolume.GetWireAngleU(), driftVolume.GetWireAngleV(), driftVolume.GetSigmaUVZ())));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, configFileName));
    }
}

} // namespace lar_pandora
