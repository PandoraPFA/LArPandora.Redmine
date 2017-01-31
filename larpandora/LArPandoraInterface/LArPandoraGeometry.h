/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.h
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#ifndef LAR_PANDORA_GEOMETRY_H
#define LAR_PANDORA_GEOMETRY_H 1

#include <string>

namespace pandora { class TiXmlElement; }

namespace lar_pandora
{

/**
 *  @brief  drift volume class to hold properties of drift volume
 */
class LArTpcVolume
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  cryostat  the cryostat ID
     *  @param  tpc       the tpc ID
     */
    LArTpcVolume(const unsigned int cryostat, unsigned int tpc);

    /**
     *  @brief  return cryostat ID
     */
    unsigned int GetCryostat() const;

    /**
     *  @brief  return tpc ID
     */
    unsigned int GetTpc() const;

private:
    unsigned int  m_cryostat;
    unsigned int  m_tpc;
};

typedef std::vector<LArTpcVolume> LArTpcVolumeList;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  drift volume class to hold properties of drift volume
 */
class LArDriftVolume
{
public:
    /**
     *  @brief  Constructor
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
     *  @param  tpcVolumeList    input list of TPC volumes
     */
    LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
        const double wirePitchU, const double wirePitchV, const double wirePitchW, const double wireAngleU, const double wireAngleV, 
        const double centerX, const double centerY, const double centerZ, const double widthX, const double widthY, const double widthZ,
        const double sigmaUVZ, const LArTpcVolumeList tpcVolumeList);

    /**
     *  @brief Return unqiue ID
     */
    unsigned int GetVolumeID() const;  

    /**
     *  @brief Return drift direction (true if positive)
     */
    bool IsPositiveDrift() const;

    /**
     *  @brief Return wire pitch in U view
     */
    double GetWirePitchU() const;

    /**
     *  @brief Return wire pictch in V view
     */
    double GetWirePitchV() const;

    /**
     *  @brief Return wire pitch in W view
     */
    double GetWirePitchW() const;

    /**
     *  @brief Return wire angle in U view (Pandora convention)
     */
    double GetWireAngleU() const;

    /**
     *  @brief Return wire angle in V view (Pandora convention)
     */
    double GetWireAngleV() const;

    /**
     *  @brief Return X position at centre of drift volume
     */
    double GetCenterX() const;

    /**
     *  @brief Return Y position at centre of drift volume
     */
    double GetCenterY() const;

    /**
     *  @brief Return Z position at centre of drift volume
     */
    double GetCenterZ() const;

    /**
     *  @brief Return X span of drift volume
     */
    double GetWidthX() const;

    /**
     *  @brief Return Y span of drift volume
     */
    double GetWidthY() const;

    /**
     *  @brief Return Z span of drift volume
     */
    double GetWidthZ() const;

    /**
     *  @brief Return sigmaUVZ parameter (used for matching views)
     */
    double GetSigmaUVZ() const;

    /**
     *  @brief Return list of tpc volumes associated with this drift volume
     */
    const LArTpcVolumeList &GetTpcVolumeList() const;

private:

    unsigned int  m_volumeID;
    bool          m_isPositiveDrift;
    double        m_wirePitchU;
    double        m_wirePitchV;
    double        m_wirePitchW;
    double        m_wireAngleU;
    double        m_wireAngleV;
    double        m_centerX;
    double        m_centerY;
    double        m_centerZ;
    double        m_widthX;
    double        m_widthY;
    double        m_widthZ;
    double        m_sigmaUVZ;

    LArTpcVolumeList m_tpcVolumeList;
};

typedef std::vector<LArDriftVolume> LArDriftVolumeList;
typedef std::map<unsigned int, LArDriftVolume> LArDriftVolumeMap;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArPandoraGeometry class
 */
class LArPandoraGeometry
{
public:
    /**
     *  @brief  This method will group TPCs into "drift volumes" (these are regions of the detector that share a common drift direction, 
     *          common range of X coordinates, and common detector parameters such as wire pitch and wire angle).
     * 
     *  @param  driftVolumeList to receive the populated drift volume list
     */
    static void LoadGeometry(LArDriftVolumeList &driftVolumeList);
 
    /**
     *  @brief  Print out the list of drift volumes
     * 
     *  @param  driftVolumeList the drift volume list
     */
    static void PrintGeometry(const LArDriftVolumeList &driftVolumeList);

    /**
     *  @brief  Write out the list of drift volumes to xml file
     * 
     *  @param  xmlFileName the xml file name
     *  @param  driftVolumeList the drift volume list
     */
    static void WriteGeometry(const std::string &xmlFileName, const LArDriftVolumeList &driftVolumeList);

private:
    /**
     *  @brief  Write an element (converted to a string) under a parent xml element
     * 
     *  @param  pParentElement the parent xml element
     *  @param  elementName the element xml name
     *  @param  elementValue the element value, converted to a string
     */
    static void WriteElement(pandora::TiXmlElement *const pParentElement, const std::string &elementName, const std::string &elementValue);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
 
inline LArTpcVolume::LArTpcVolume(const unsigned int cryostat, unsigned int tpc) :
    m_cryostat(cryostat), m_tpc(tpc)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArTpcVolume::GetCryostat() const
{
    return m_cryostat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArTpcVolume::GetTpc() const
{
    return m_tpc;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArDriftVolume::GetVolumeID() const
{ 
    return m_volumeID; 
}
   
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool LArDriftVolume::IsPositiveDrift() const
{
    return m_isPositiveDrift;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------      

inline double LArDriftVolume::GetWirePitchU() const
{ 
    return m_wirePitchU; 
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWirePitchV() const
{ 
    return m_wirePitchV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWirePitchW() const
{
    return m_wirePitchW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWireAngleU() const
{
    return m_wireAngleU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWireAngleV() const
{
    return m_wireAngleV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetCenterX() const
{
    return m_centerX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetCenterY() const
{ 
    return m_centerY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetCenterZ() const
{
    return m_centerZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWidthX() const
{
    return m_widthX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWidthY() const
{ 
    return m_widthY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetWidthZ() const
{
    return m_widthZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDriftVolume::GetSigmaUVZ() const
{
    return m_sigmaUVZ;
}

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_GEOMETRY_H
