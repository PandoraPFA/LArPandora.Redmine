/**
 *  @file   larpandora/LArPandoraInterface/LArPandoraGeometry.h
 *
 *  @brief  Helper functions for extracting detector geometry for use in reconsruction
 */

#ifndef LAR_PANDORA_GEOMETRY_H
#define LAR_PANDORA_GEOMETRY_H 1

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
    LArDetectorGap(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2);

    /**
     *  @brief Get lower X coordinate
     */
    float GetX1() const;

    /**
     *  @brief Get lower y coordinate
     */
    float GetY1() const;

    /**
     *  @brief Get lower Z coordinate
     */
    float GetZ1() const;

    /**
     *  @brief Get upper X coordinate
     */
    float GetX2() const;

    /**
     *  @brief Get upper Y coordinate
     */
    float GetY2() const;

    /**
     *  @brief Get upper Z coordinate
     */
    float GetZ2() const;

private:
    float   m_x1;
    float   m_y1;
    float   m_z1;
    float   m_x2;
    float   m_y2;
    float   m_z2;
};

typedef std::vector<LArDetectorGap> LArDetectorGapList;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  daughter drift volume class to hold properties of daughter drift volumes
 */
class LArDaughterDriftVolume
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  cryostat  the cryostat ID
     *  @param  tpc       the tpc ID
     */
    LArDaughterDriftVolume(const unsigned int cryostat, const unsigned int tpc);

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

typedef std::vector<LArDaughterDriftVolume> LArDaughterDriftVolumeList;

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
     *  @param  thetaU           wire angle to vertical (U)
     *  @param  thetaV           wire angle to vertical (V)
     *  @param  sigmaUVZ         matching between views
     *  @param  tpcVolumeList    input list of TPC volumes
     */
    LArDriftVolume(const unsigned int volumeID, const bool isPositiveDrift,
        const float wirePitchU, const float wirePitchV, const float wirePitchW, const float wireAngleU, const float wireAngleV,
        const float centerX, const float centerY, const float centerZ, const float widthX, const float widthY, const float widthZ,
        const float sigmaUVZ, const LArDaughterDriftVolumeList &tpcVolumeList);

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
    float GetWirePitchU() const;

    /**
     *  @brief Return wire pictch in V view
     */
    float GetWirePitchV() const;

    /**
     *  @brief Return wire pitch in W view
     */
    float GetWirePitchW() const;

    /**
     *  @brief Return wire angle in U view (Pandora convention)
     */
    float GetWireAngleU() const;

    /**
     *  @brief Return wire angle in V view (Pandora convention)
     */
    float GetWireAngleV() const;

    /**
     *  @brief Return X position at centre of drift volume
     */
    float GetCenterX() const;

    /**
     *  @brief Return Y position at centre of drift volume
     */
    float GetCenterY() const;

    /**
     *  @brief Return Z position at centre of drift volume
     */
    float GetCenterZ() const;

    /**
     *  @brief Return X span of drift volume
     */
    float GetWidthX() const;

    /**
     *  @brief Return Y span of drift volume
     */
    float GetWidthY() const;

    /**
     *  @brief Return Z span of drift volume
     */
    float GetWidthZ() const;

    /**
     *  @brief Return sigmaUVZ parameter (used for matching views)
     */
    float GetSigmaUVZ() const;

    /**
     *  @brief Return list of daughter drift volumes associated with this drift volume
     */
    const LArDaughterDriftVolumeList &GetTpcVolumeList() const;

private:
    unsigned int    m_volumeID;
    bool            m_isPositiveDrift;
    float           m_wirePitchU;
    float           m_wirePitchV;
    float           m_wirePitchW;
    float           m_wireAngleU;
    float           m_wireAngleV;
    float           m_centerX;
    float           m_centerY;
    float           m_centerZ;
    float           m_widthX;
    float           m_widthY;
    float           m_widthZ;
    float           m_sigmaUVZ;

    LArDaughterDriftVolumeList m_tpcVolumeList;
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
    static void LoadGeometry(const Settings &settings, LArDriftVolumeList &driftVolumeList);

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
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArDetectorGap::LArDetectorGap(const float x1, const float y1, const float z1, const float x2, const float y2, const float z2) :
    m_x1(x1), m_y1(y1), m_z1(z1), m_x2(x2), m_y2(y2), m_z2(z2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDetectorGap::GetX1() const
{
    return m_x1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDetectorGap::GetY1() const
{
    return m_y1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDetectorGap::GetZ1() const
{
    return m_z1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDetectorGap::GetX2() const
{
    return m_x2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDetectorGap::GetY2() const
{
    return m_y2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDetectorGap::GetZ2() const
{
    return m_z2;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArDaughterDriftVolume::LArDaughterDriftVolume(const unsigned int cryostat, const unsigned int tpc) :
    m_cryostat(cryostat), m_tpc(tpc)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArDaughterDriftVolume::GetCryostat() const
{
    return m_cryostat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArDaughterDriftVolume::GetTpc() const
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

inline float LArDriftVolume::GetWirePitchU() const
{
    return m_wirePitchU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWirePitchV() const
{
    return m_wirePitchV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWirePitchW() const
{
    return m_wirePitchW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWireAngleU() const
{
    return m_wireAngleU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWireAngleV() const
{
    return m_wireAngleV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetCenterX() const
{
    return m_centerX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetCenterY() const
{
    return m_centerY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetCenterZ() const
{
    return m_centerZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWidthX() const
{
    return m_widthX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWidthY() const
{
    return m_widthY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetWidthZ() const
{
    return m_widthZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArDriftVolume::GetSigmaUVZ() const
{
    return m_sigmaUVZ;
}

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_GEOMETRY_H
