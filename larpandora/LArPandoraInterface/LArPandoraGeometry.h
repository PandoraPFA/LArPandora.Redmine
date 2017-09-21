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
class LArDetectorGap
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  x1 lower X coordinate
     *  @param  y1 lower Y coordinate
     *  @param  z1 lower Z coordinate
     *  @param  x2 upper X coordinate
     *  @param  y2 upper Y coordinate
     *  @param  z2 upper Z coordinate
     */
    LArDetectorGap(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2);

    /**
     *  @brief Get lower X coordinate
     */
    double GetX1() const;

    /**
     *  @brief Get lower y coordinate
     */
    double GetY1() const;

    /**
     *  @brief Get lower Z coordinate
     */
    double GetZ1() const;

    /**
     *  @brief Get upper X coordinate
     */
    double GetX2() const;

    /**
     *  @brief Get upper Y coordinate
     */
    double GetY2() const;

    /**
     *  @brief Get upper Z coordinate
     */
    double GetZ2() const;

private:
    double  m_x1;
    double  m_y1;
    double  m_z1;
    double  m_x2;
    double  m_y2;
    double  m_z2;
};

typedef std::vector<LArDetectorGap> LArDetectorGapList;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

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
    LArTpcVolume(const unsigned int cryostat, const unsigned int tpc);

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
     *  @brief Return unique ID
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
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        bool            m_globalCoordinates;        ///< Use global coordinate system
        bool            m_globalDriftVolume;        ///< Use global drift volume
        bool            m_printGeometry;            ///< Print geometry to screen
        std::string     m_geometryXmlFileName;      ///< File to which to write compact geometry xml
    };

    /**
     *  @brief Load the 2D gaps that go with the chosen geometry
     *
     *  @param settings the input configuration
     *  @param listOfGaps the output list of 2D gaps.
     */
    static void LoadDetectorGaps(const Settings &settings, LArDetectorGapList &listOfGaps);

    /**
     *  @brief Load drift volume geometry
     *
     *  @param settings the input configuration
     *  @param driftVolumeList the output list of drift volumes
     *  @param driftVolumeMap the output mapping between cryostat/tpc and drift volumes
     */
    static void LoadGeometry(const Settings &settings, LArDriftVolumeList &driftVolumeList, LArDriftVolumeMap &driftVolumeMap);

    /**
     *  @brief  Print out the list of drift volumes
     *
     *  @param  settings the input configuration (so we can print out this as well)
     *  @param  driftVolumeList the drift volume list
     */
    static void PrintGeometry(const Settings &settings, const LArDriftVolumeList &driftVolumeList);

    /**
     *  @brief  Write out the list of drift volumes to xml file
     *
     *  @param  xmlFileName the xml file name
     *  @param  driftVolumeList the drift volume list
     */
    static void WriteGeometry(const std::string &xmlFileName, const LArDriftVolumeList &driftVolumeList);

    /**
     *  @brief  Get drift volume ID from a specified cryostat/tpc pair
     *
     *  @param  driftVolumeMap the output mapping between cryostat/tpc and drift volumes
     *  @param  cstat the input cryostat unique ID
     *  @param  tpc the input tpc unique ID
     */
    static unsigned int GetVolumeID(const LArDriftVolumeMap &driftVolumeMap, const unsigned int cstat, const unsigned int tpc);

    /**
     *  @brief  Convert to global coordinate system
     *
     *  @param  cstat the input cryostat
     *  @param  tpc the input tpc
     *  @param  hit_View the input view
     */
    static geo::View_t GetGlobalView(const unsigned int cstat, const unsigned int tpc, const geo::View_t hit_View);

private:

    /**
     *  @brief  Generate a unique identifier for each TPC
     *
     *  @param  cstat the input cryostat
     *  @param  tpc the input tpc
     */
    static unsigned int GetTpcID(const unsigned int cstat, const unsigned int tpc);

    /**
     *  @brief  Return whether U/V should be switched in global coordinate system for this cryostat/tpc
     *
     *  @param  cstat the input cryostat
     *  @param  tpc the input tpc
     */
    static bool ShouldSwitchUV(const unsigned int cstat, const unsigned int tpc);

    /**
     *  @brief  Return whether U/V should be switched in global coordinate system for this drift direction
     *
     *  @param  isPositiveDrift the drift direction
     */
    static bool ShouldSwitchUV(const bool isPositiveDrift);

    /**
     *  @brief  This method will group TPCs into drift volumes (these are regions of the detector that share a common drift direction,
     *          common range of X coordinates, and common detector parameters such as wire pitch and wire angle).
     *
     *  @param  driftVolumeList to receive the populated drift volume list
     */
    static void LoadGeometry(LArDriftVolumeList &driftVolumeList);

    /**
     *  @brief  This method will create one or more global volumes (these are conglomeration of drift volumes that share a common drift
     *          orientation along the X-axis, have parallel or near-parallel wire angles, and similar wire pitches)
     *
     *  @param  driftVolumeList to receive the input drift volume list
     *  @param  parentVolumeList to receive the output parent drift volume list
     */
    static void LoadGlobalParentGeometry(const LArDriftVolumeList &driftVolumeList, LArDriftVolumeList &parentVolumeList);

    /**
     *  @brief  This method will create one or more daughter volumes (these share a common drift orientation along the X-axis,
     *          have parallel or near-parallel wire angles, and similar wire pitches)
     *
     *  @param  driftVolumeList to receive the input drift volume list
     *  @param  parentVolumeList to receive the output daughter drift volume list
     */
    static void LoadGlobalDaughterGeometry(const LArDriftVolumeList &driftVolumeList, LArDriftVolumeList &daughterVolumeList);
    /**
     *  @brief  Write a double element, above standard precision (required for e.g. wire angle values), under a parent xml element
     *
     *  @param  pParentElement the parent xml element
     *  @param  elementName the element xml name
     *  @param  value the element value
     */
    static void WritePrecisionElement(pandora::TiXmlElement *const pParentElement, const std::string &elementName, const double value);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArDetectorGap::LArDetectorGap(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2) :
    m_x1(x1), m_y1(y1), m_z1(z1), m_x2(x2), m_y2(y2), m_z2(z2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDetectorGap::GetX1() const
{
    return m_x1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDetectorGap::GetY1() const
{
    return m_y1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDetectorGap::GetZ1() const
{
    return m_z1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDetectorGap::GetX2() const
{
    return m_x2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDetectorGap::GetY2() const
{
    return m_y2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArDetectorGap::GetZ2() const
{
    return m_z2;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTpcVolume::LArTpcVolume(const unsigned int cryostat, const unsigned int tpc) :
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
