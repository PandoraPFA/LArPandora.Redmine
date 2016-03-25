/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraModule_module.cc
 *
 *  @brief  LArPandora producer module for MicroBooNE
 *
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "LArStitching/MultiPandoraApi.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"

#include <string>

namespace lar_pandora
{

/**
 *  @brief  LArPandoraModule class
 */
class LArPandoraModule : public LArPandora
{
public: 

    /**
     *  @brief  drift volume class to hold properties of drift volume
     */
    class LArDriftVolume
    {
      public:
        /**
         *  @brief  Default constructor
         *
         *  @param  volumeID         unique ID number
         *  @param  isPositiveDrift  direction of drift 
         *  @param  wirePitchU       wire pitch (U view)
         *  @param  wirePitchV       wire pitch (V view)
         *  @param  wirePitchW       wire pitch (W view)
         *  @param  wireAngleU       wire angle (U view)
         *  @param  wireAngleV       wire angle (V view)
         *  @param  centerX          centre of volume (X)
         *  @param  centerY          centre of volume (Y)
         *  @param  centerZ          centre of volume (Z)
         *  @param  widthX           width of volume (X)
         *  @param  widthY           width of volume (Y)
         *  @param  widthZ           width of volume (Z)
         *  @param  sigmaUVZ         matching between views
         */
        LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
            const double wirePitchU, const double wirePitchV, const double wirePitchW, const double wireAngleU, const double wireAngleV, 
	    const double centerX, const double centerY, const double centerZ, const double widthX, const double widthY, const double widthZ,
            const double sigmaUVZ);

        const unsigned int  m_volumeID;
        const bool          m_isPositiveDrift;
        const double        m_wirePitchU;
        const double        m_wirePitchV; 
        const double        m_wirePitchW;
        const double        m_wireAngleU;
        const double        m_wireAngleV;
        const double        m_centerX;
        const double        m_centerY;
        const double        m_centerZ;
        const double        m_widthX;
        const double        m_widthY;
        const double        m_widthZ;
        const double        m_sigmaUVZ;
    };

    typedef std::vector<LArDriftVolume> LArDriftVolumeList;
    typedef std::map<unsigned int, LArDriftVolume> LArDriftVolumeMap;
    

    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    LArPandoraModule(fhicl::ParameterSet const &pset);

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
     *  @brief  method to extract geometry information from LArSoft and fill vector of drift volume objects
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

    LArDriftVolumeList  m_driftVolumeList;
    LArDriftVolumeMap   m_driftVolumeMap;
 
    bool     m_printDebug;       // Print manual debug messages
    double   m_maxDeltaTheta;    // Allowed variation in common wire angle within a drift volume

    bool     m_useShortVolume;   // Historical DUNE 35t config parameter - use short drift volume
    bool     m_useLongVolume;    // Historical DUNE 35t config parameter - use long drift volume

    bool     m_useLeftVolume;    // Historical DUNE 4-APA config parameter - use left drift volume
    bool     m_useRightVolume;   // Historical DUNE 4-APA config parameter - use right drift volume
};

DEFINE_ART_MODULE(LArPandoraModule)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib/exception.h"

#include "larcore/Geometry/Geometry.h"

#include "Api/PandoraApi.h"

#include "LArContent.h"
#include "LArPlugins/LArPseudoLayerPlugin.h"
#include "LArPlugins/LArRotationalTransformationPlugin.h"

namespace lar_pandora
{

LArPandoraModule::LArPandoraModule(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
    m_printDebug = false;
    m_maxDeltaTheta = 0.01;

    m_useShortVolume = pset.get<bool>("UseShortVolume", true);
    m_useLongVolume = pset.get<bool>("UseLongVolume", true);

    m_useLeftVolume = pset.get<bool>("UseLeftVolume", true);
    m_useRightVolume = pset.get<bool>("UseRightVolume", true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArPandoraModule::GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const
{
    LArDriftVolumeMap::const_iterator iter = m_driftVolumeMap.find(this->GetTpcIdNumber(cryostat, tpc));

    if (m_driftVolumeMap.end() == iter)
        throw cet::exception("LArPandora") << " Throwing exception - found a TPC that doesn't belong to a drift volume";
 
    const LArDriftVolume &driftVolume = iter->second;
    return driftVolume.m_volumeID;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArPandoraModule::GetTpcIdNumber(const unsigned int cryostat, const unsigned int tpc) const
{
    return 10000 * cryostat + tpc;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraModule::LoadGeometry()
{ 
    // This method will group TPCs into "drift volumes" (these are regions of the detector that share a common drift direction, 
    // common range of X coordinates, and common detector parameters such as wire pitch and wire angle).
   
    mf::LogDebug("LArPandora") << " *** LArPandoraModule::LoadGeometry() *** " << std::endl;

    if (!m_driftVolumeList.empty() || !m_driftVolumeMap.empty())
        throw cet::exception("LArPandora") << " Throwing exception - detector geoemetry has already been loaded ";

    typedef std::set<unsigned int> UIntSet;

    // Loading Geometry Service
    art::ServiceHandle<geo::Geometry> theGeometry;

    mf::LogDebug("LArPandora") << " Using Geometry: " << theGeometry->DetectorName() << std::endl;

    const unsigned int wirePlanes(theGeometry->MaxPlanes());

    const double wirePitchU(theGeometry->WirePitch(geo::kU));
    const double wirePitchV(theGeometry->WirePitch(geo::kV));
    const double wirePitchW((wirePlanes > 2) ? theGeometry->WirePitch(geo::kW) : 0.5 * (wirePitchU + wirePitchV));

    // Loop over cryostats
    for (unsigned int icstat = 0; icstat < theGeometry->Ncryostats(); ++icstat)
    {
        UIntSet cstatList;

        // First loop over TPCs
        for (unsigned int itpc1 = 0; itpc1 < theGeometry->NTPC(icstat); ++itpc1)
        {      
	    if (cstatList.end() != cstatList.find(itpc1))
	        continue;

            // Use this TPC to seed a drift volume
            const geo::TPCGeo &theTpc1(theGeometry->TPC(itpc1, icstat));
            cstatList.insert(itpc1);

            const double wireAngleU(0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kU, itpc1, icstat));
	    const double wireAngleV((0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kV, itpc1, icstat)) * -1.f);
            const double wireAngleZ((wirePlanes > 2) ? (0.5f * M_PI - theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat)) : 0.0);

            if (std::fabs(wireAngleZ) > m_maxDeltaTheta)
                throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

            double localCoord1[3] = {0.,0.,0.};
            double worldCoord1[3] = {0.,0.,0.};
            theTpc1.LocalToWorld(localCoord1, worldCoord1);

            const double min1(worldCoord1[0] - 0.5 * theTpc1.ActiveHalfWidth());
            const double max1(worldCoord1[0] + 0.5 * theTpc1.ActiveHalfWidth());

            double driftMinX(worldCoord1[0] - theTpc1.ActiveHalfWidth());
            double driftMaxX(worldCoord1[0] + theTpc1.ActiveHalfWidth());
            double driftMinY(worldCoord1[1] - theTpc1.ActiveHalfHeight());
            double driftMaxY(worldCoord1[1] + theTpc1.ActiveHalfHeight());
            double driftMinZ(worldCoord1[2] - 0.5 * theTpc1.ActiveLength());
            double driftMaxZ(worldCoord1[2] + 0.5 * theTpc1.ActiveLength());

            const bool isPositiveDrift(theTpc1.DriftDirection() == geo::kPosX);

            UIntSet tpcList;
            tpcList.insert(itpc1);

            // Now identify the other TPCs associated with this drift volume
            for (unsigned int itpc2 = itpc1+1; itpc2 < theGeometry->NTPC(icstat); ++itpc2)
            {
                if (cstatList.end() != cstatList.find(itpc2))
	            continue;
         
                const geo::TPCGeo &theTpc2(theGeometry->TPC(itpc2, icstat));

                if (theTpc1.DriftDirection() != theTpc2.DriftDirection())
		    continue;

                const double dThetaU(theGeometry->WireAngleToVertical(geo::kU, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kU, itpc2, icstat));
                const double dThetaV(theGeometry->WireAngleToVertical(geo::kV, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kV, itpc2, icstat));
                const double dThetaW((wirePlanes > 2) ? (theGeometry->WireAngleToVertical(geo::kW, itpc1, icstat) - theGeometry->WireAngleToVertical(geo::kW, itpc2, icstat)) : 0.0);

                if (dThetaU > m_maxDeltaTheta || dThetaV > m_maxDeltaTheta || dThetaW > m_maxDeltaTheta)
		   continue;

                double localCoord2[3] = {0.,0.,0.};
                double worldCoord2[3] = {0.,0.,0.};
                theTpc2.LocalToWorld(localCoord2, worldCoord2);
                const double min2(worldCoord2[0] - 0.5 * theTpc2.ActiveHalfWidth());
                const double max2(worldCoord2[0] + 0.5 * theTpc2.ActiveHalfWidth());

                if ((min2 > max1) || (min1 > max2))
		    continue;

                cstatList.insert(itpc2);
                tpcList.insert(itpc2);

                driftMinX = std::min(driftMinX, worldCoord2[0] - theTpc2.ActiveHalfWidth());
                driftMaxX = std::max(driftMaxX, worldCoord2[0] + theTpc2.ActiveHalfWidth());
                driftMinY = std::min(driftMinY, worldCoord2[1] - theTpc2.ActiveHalfHeight());
                driftMaxY = std::max(driftMaxY, worldCoord2[1] + theTpc2.ActiveHalfHeight());
                driftMinZ = std::min(driftMinZ, worldCoord2[2] - 0.5 * theTpc2.ActiveLength());
                driftMaxZ = std::max(driftMaxZ, worldCoord2[2] + 0.5 * theTpc2.ActiveLength());
	    }

            // create new drift volume
            const LArDriftVolume driftVolume(m_driftVolumeList.size(), isPositiveDrift,
                wirePitchU, wirePitchV, wirePitchW, wireAngleU, wireAngleV,
	        0.5 * (driftMaxX + driftMinX), 0.5 * (driftMaxY + driftMinY), 0.5 * (driftMaxZ + driftMinZ),
		(driftMaxX - driftMinX), (driftMaxY - driftMinY), (driftMaxZ - driftMinZ),
		(wirePitchU + wirePitchV + wirePitchW + 0.1));

            m_driftVolumeList.push_back(driftVolume);
            
            if (m_printDebug)
	    {
	        std::cout << " *** New Drift Volume *** " << std::endl;
	        std::cout << "  ID = " << driftVolume.m_volumeID << std::endl;
                std::cout << "  isPositiveDrift = " <<driftVolume.m_isPositiveDrift << std::endl;
                std::cout << "  m_wirePitchU = " << driftVolume.m_wirePitchU << std::endl;
                std::cout << "  m_wirePitchV = " << driftVolume.m_wirePitchV << std::endl;
                std::cout << "  m_wirePitchW = " << driftVolume.m_wirePitchW << std::endl;
                std::cout << "  m_wireAngleU = " << driftVolume.m_wireAngleU << std::endl;
                std::cout << "  m_wireAngleV = " << driftVolume.m_wireAngleV << std::endl;
	        std::cout << "  m_centerX = " << driftVolume.m_centerX << std::endl;
	        std::cout << "  m_centerY = " << driftVolume.m_centerY << std::endl;
	        std::cout << "  m_centerZ = " << driftVolume.m_centerZ << std::endl;
	        std::cout << "  m_widthX = " << driftVolume.m_widthX << std::endl;
	        std::cout << "  m_widthY = " << driftVolume.m_widthY << std::endl;
	        std::cout << "  m_widthZ = " << driftVolume.m_widthZ << std::endl;
		std::cout << "  m_sigmaUVZ = " << driftVolume.m_sigmaUVZ << std::endl;
		std::cout << "  TPC List: " << std::endl;
            }

            for(UIntSet::const_iterator iter = tpcList.begin(), iterEnd = tpcList.end(); iter != iterEnd; ++iter)
	    {
                m_driftVolumeMap.insert(LArDriftVolumeMap::value_type(this->GetTpcIdNumber(icstat, *iter), driftVolume));

                if (m_printDebug)
		    std::cout << "   [" << icstat << "][" << *iter << "]" << std::endl;
	    }
	}
    }

    if (m_driftVolumeList.empty() || m_driftVolumeMap.empty())
        throw cet::exception("LArPandora") << " Failed to find any drift volumes in this detector geometry ";
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraModule::CreatePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** LArPandoraModule::CreatePandoraInstances(...) [BEGIN] *** " << std::endl;
    
    // Load geometry
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

    mf::LogDebug("LArPandora") << " *** LArPandoraModule::CreatePandoraInstances(...) [END] *** " << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraModule::CreatePrimaryPandoraInstance(const std::string &configFileName)
{
    mf::LogDebug("LArPandora") << " *** LArPandoraModule::CreatePrimaryPandoraInstance(...) *** " << std::endl;

    const LArDriftVolume &driftVol = *(m_driftVolumeList.begin());

    if (m_driftVolumeList.size() > 1)
    {
        m_pPrimaryPandora = this->CreateNewPandora();
        MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, configFileName));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora,
	    new lar_content::LArPseudoLayerPlugin(driftVol.m_wirePitchU, driftVol.m_wirePitchV, driftVol.m_wirePitchW)));
    }

    else
    {
        m_pPrimaryPandora = this->CreateNewPandora();
        MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
        MultiPandoraApi::SetVolumeInfo(m_pPrimaryPandora, new VolumeInfo(0, "driftVolume",
            pandora::CartesianVector(driftVol.m_centerX, driftVol.m_centerY, driftVol.m_centerZ), driftVol.m_isPositiveDrift));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora,
	    new lar_content::LArPseudoLayerPlugin(driftVol.m_wirePitchU, driftVol.m_wirePitchV, driftVol.m_wirePitchW)));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*m_pPrimaryPandora,
            new lar_content::LArRotationalTransformationPlugin(driftVol.m_wireAngleU, driftVol.m_wireAngleV, driftVol.m_sigmaUVZ)));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, configFileName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraModule::CreateDaughterPandoraInstances(const std::string &configFileName)
{
    mf::LogDebug("LArPandora") << " *** LArPandoraModule::CreateDaughterPandoraInstance(...) *** " << std::endl;

    if (!m_pPrimaryPandora)
      throw cet::exception("LArPandora") << " Throwing exception - trying to create daughter Pandora instances in absence of primary instance ";

    for (LArDriftVolumeList::const_iterator iter = m_driftVolumeList.begin(), iterEnd = m_driftVolumeList.end(); iter != iterEnd; ++iter)
    {
        const LArDriftVolume &driftVol = *iter;

        // Check historical DUNE config parameters
        if ((2 == m_driftVolumeList.size()) &&
            ((true == driftVol.m_isPositiveDrift && (false == m_useShortVolume || false == m_useRightVolume)) ||
             (false == driftVol.m_isPositiveDrift && (false == m_useLongVolume || false == m_useLeftVolume))))
	    continue;

	mf::LogDebug("LArPandora") << " Creating Pandora Daughter Instance: [" << driftVol.m_volumeID << "]" << std::endl;

        const unsigned int volumeIdNumber(driftVol.m_volumeID);
        std::ostringstream volumeIdString("driftVolume_");
        volumeIdString << volumeIdNumber;

        const pandora::Pandora *const pPandora = this->CreateNewPandora();
        MultiPandoraApi::AddDaughterPandoraInstance(m_pPrimaryPandora, pPandora);
        MultiPandoraApi::SetVolumeInfo(pPandora, new VolumeInfo(volumeIdNumber, volumeIdString.str(),
            pandora::CartesianVector(driftVol.m_centerX, driftVol.m_centerY, driftVol.m_centerZ), driftVol.m_isPositiveDrift));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora,
	    new lar_content::LArPseudoLayerPlugin(driftVol.m_wirePitchU, driftVol.m_wirePitchV, driftVol.m_wirePitchW)));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora,
	    new lar_content::LArRotationalTransformationPlugin(driftVol.m_wireAngleU, driftVol.m_wireAngleV, driftVol.m_sigmaUVZ)));
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, configFileName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArPandoraModule::LArDriftVolume::LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
  const double wirePitchU, const double wirePitchV, const double wirePitchW, const double wireAngleU, const double wireAngleV,
  const double centerX, const double centerY, const double centerZ, const double widthX, const double widthY, const double widthZ,
  const double sigmaUVZ) :
    m_volumeID(volumeID), m_isPositiveDrift(isPositiveDrift), m_wirePitchU(wirePitchU), m_wirePitchV(wirePitchV), m_wirePitchW(wirePitchW),
    m_wireAngleU(wireAngleU), m_wireAngleV(wireAngleV), m_centerX(centerX), m_centerY(centerY), m_centerZ(centerZ),
    m_widthX(widthX), m_widthY(widthY), m_widthZ(widthZ), m_sigmaUVZ(sigmaUVZ)
{
}

} // namespace lar_pandora
