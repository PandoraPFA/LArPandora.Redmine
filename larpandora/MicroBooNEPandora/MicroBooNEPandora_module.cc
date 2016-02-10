/**
 *  @file   larpandora/MicroBooNEPandora/MicroBooNEPandora_module.cc
 *
 *  @brief  LArPandora producer module for MicroBooNE
 *
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "LArStitching/MultiPandoraApi.h"

#include "LArPandoraInterface/LArPandora.h"

#include <string>

namespace lar_pandora
{

/**
 *  @brief  MicroBooNEPandora class
 */
class MicroBooNEPandora : public LArPandora
{
public: 
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    MicroBooNEPandora(fhicl::ParameterSet const &pset);

private:
    void CreatePandoraInstances();
    int GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const;

    /**
     *  @brief  Create primary pandora instance
     *
     *  @param  configFileName the pandora settings config file name
     */
    void CreatePrimaryPandoraInstance(const std::string &configFileName);
};

DEFINE_ART_MODULE(MicroBooNEPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib/exception.h"

#include "Geometry/Geometry.h"

#include "Api/PandoraApi.h"

#include "LArContent.h"

#include "MicroBooNEPandora/MicroBooNEPseudoLayerPlugin.h"
#include "MicroBooNEPandora/MicroBooNETransformationPlugin.h"

namespace lar_pandora
{

MicroBooNEPandora::MicroBooNEPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

int MicroBooNEPandora::GetVolumeIdNumber(const unsigned int /*cryostat*/, const unsigned int /*tpc*/) const
{
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::CreatePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** MicroBooNEPandora::CreatePandoraInstances(...) *** " << std::endl;

    art::ServiceHandle<geo::Geometry> theGeometry;
    if (std::string::npos == theGeometry->DetectorName().find("microboone"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " MicroBooNEPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    std::string configFileName;
    cet::search_path sp("FW_SEARCH_PATH");
    if (!sp.find_file(m_configFile, configFileName))
    {
        mf::LogError("LArPandora") << "   Failed to find: " << m_configFile << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);              
    }                                                                                    

    this->CreatePrimaryPandoraInstance(configFileName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::CreatePrimaryPandoraInstance(const std::string &configFileName)
{
    m_pPrimaryPandora = this->CreateNewPandora();
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
    MultiPandoraApi::SetVolumeInfo(m_pPrimaryPandora, new VolumeInfo(0, "microboone", pandora::CartesianVector(0.f, 0.f, 0.f), true));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora, new lar_pandora::MicroBooNEPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*m_pPrimaryPandora, new lar_pandora::MicroBooNETransformationPlugin));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, configFileName));
}

} // namespace lar_pandora
