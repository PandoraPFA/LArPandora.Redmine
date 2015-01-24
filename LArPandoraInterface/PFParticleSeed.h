/**
 *  @file   larpandora/LArPandoraInterface/PFParticleSeed.h
 *
 *  @brief  Class to hold start/end position and directions of PFParticles
 *
 */
#ifndef PFPARTICLE_SEED_H
#define PFPARTICLE_SEED_H 1

// Local includes
#include "LArPandoraCollector.h"

// Pandora includes
#include "Helpers/ClusterFitHelper.h"

// ROOT includes
#include "TVector3.h"

namespace lar_pandora
{

/**
 *  @brief  PFParticleTrajectoryPoint class
 */
class PFParticleTrajectoryPoint
{
  public:  
    /**
     *  @brief  Constructor
     *
     *  @param  position
     *  @param  direction
     *  @param  displacement
     */
    PFParticleTrajectoryPoint(const pandora::CartesianVector &position, const pandora::CartesianVector &direction, 
        const float displacement);

    pandora::CartesianVector  m_position; 
    pandora::CartesianVector  m_direction; 
    float                     m_displacement; 

    bool operator < (const PFParticleTrajectoryPoint& rhs) const
    {
        return (m_displacement < rhs.m_displacement);
    }
};

typedef std::vector<PFParticleTrajectoryPoint> PFParticleTrajectoryPointList;


/**
 *  @brief  PFParticleSeed class
 */
class PFParticleSeed
{
public:

    /**
     *  @brief  Constructor
     *
     *  @param  SpacePointVector
     */
    PFParticleSeed(const SpacePointVector &spacepoints);

    /**
     *  @brief  Constructor
     *
     *  @param  SpacePointVector
     */
    PFParticleSeed(const pandora::CartesianPointList &pointList);

    /**
     *  @brief  Destructor
     */
    ~PFParticleSeed();

    /**
     *  @brief Initialise Particle Seed
     */
    void Initialize(const pandora::CartesianPointList &pointList);

    /**
     *  @brief Check initialisation
     */
    bool IsInitialized() const;

    /**
     *  @brief  Return vertex position 
     */
    pandora::CartesianVector GetInnerPosition() const;

    /**
     *  @brief  Return vertex direction 
     */
    pandora::CartesianVector GetInnerDirection() const;

    /**
     *  @brief  Return end position
     */
    pandora::CartesianVector GetOuterPosition()  const;

    /**
     *  @brief  Return end direction
     */
    pandora::CartesianVector GetOuterDirection() const;

private:

    bool                           m_isInitialized;    ///<
    pandora::CartesianVector       m_innerPosition;    ///<
    pandora::CartesianVector       m_innerDirection;   ///<
    pandora::CartesianVector       m_outerPosition;    ///<
    pandora::CartesianVector       m_outerDirection;   ///<
};


/**
 *  @brief  PFParticleFitter class
 */
class PFParticleFitter
{
public:

    /**
     *  @brief Build a vector of PANDORA space points from a vector of ART space points
     *
     *  @param  spacepoints  the input list of ART space points
     *  @param  pointList  the output list of PANDORA space points
     */
    static void BuildPointList(const SpacePointVector &spacepoints, pandora::CartesianPointList &pointList);

    /**
     *  @brief Build a vector of trajectory points from a vector of space points
     *
     *  @param  spacepoints  the input list of space points
     *  @param  trajectoryPointList  the output list of trajectory points
     *  @param  isCosmic  choice of vertex (true for highest vertex in Y, false for lowest vertex in Z) 
     */
    static void BuildTrajectoryPointList(const SpacePointVector &spacepoints, PFParticleTrajectoryPointList &trajectoryPointList,
        const bool isCosmic = true);

    /**
     *  @brief Find the extremal coordinates from a given set of space points
     *
     *  @param  pointList  the input list of space points
     *  @param  innerCoordinate  the inner extremal coordinate (minimum z value)
     *  @param  outerCoordinate  the outer extremal coordinate (maximum z value)
     */
    static void GetExtremalCoordinates(const pandora::CartesianPointList &pointList, pandora::CartesianVector &innerCoordinate, 
        pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief Fit the start/end region of a given set of space points
     *
     *  @param  pointList  the input list of space points
     *  @param  vertexPosition  the input seed vertex position
     *  @param  vertexDirection  the input seed vertex direction 
     *  @param  fittedDirection  the output fitted vertex direction
     */
    static void FitPoints(const pandora::CartesianPointList &pointList, const pandora::CartesianVector &vertexPosition,
        const pandora::CartesianVector &vertexDirection, pandora::CartesianVector &fittedDirection);
};

} // namespace lar_pandora

#endif // #ifndef PFPARTICLE_SEED_H
