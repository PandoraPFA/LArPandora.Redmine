/**
 *  @file   larpandora/LArPandoraInterface/MicroBooNEPandora_module.cc
 *
 *  @brief  Producer module for MicroBooNE
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local LArPandora includes
#include "larpandora/LArPandoraInterface/LArPandoraParticleCreator.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  MicroBooNEPandora class
 */
class MicroBooNEPandora : public LArPandoraParticleCreator
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    MicroBooNEPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~MicroBooNEPandora();

private:
    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;
};

DEFINE_ART_MODULE(MicroBooNEPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// Framework includes
#include "cetlib/exception.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"

// LArContent includes
#include "LArContent.h"

// Local includes
#include "larpandora/LArPandoraInterface/MicroBooNEPseudoLayerPlugin.h"
#include "larpandora/LArPandoraInterface/MicroBooNETransformationPlugin.h"

namespace lar_pandora {

MicroBooNEPandora::MicroBooNEPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{  
}

//------------------------------------------------------------------------------------------------------------------------------------------

MicroBooNEPandora::~MicroBooNEPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MicroBooNEPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** MicroBooNEPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("microboone"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " MicroBooNEPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const pandora::Pandora *pPandora = pIter->second;
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new MicroBooNEPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new MicroBooNETransformationPlugin));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int MicroBooNEPandora::GetPandoraVolumeID(const unsigned int, const unsigned int) const
{
    return 0;
}

} // namespace lar_pandora
