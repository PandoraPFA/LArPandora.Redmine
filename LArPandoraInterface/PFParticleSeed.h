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

namespace lar_pandora
{

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
     *  @brief  Destructor
     */
    ~PFParticleSeed();

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

    /**
     *  @brief Build a list of PANDORA space points from a vector of ART space points
     *
     *  @param  spacepoints  the input list of ART space points
     *  @param  pointList  the output list of PANDORA space points
     */
    void BuildPointList(const SpacePointVector &spacepoints, pandora::CartesianPointList &pointList) const;

    /**
     *  @brief Find the extremal coordinates from a given set of space points
     *
     *  @param  pointList  the input list of space points
     *  @param  innerCoordinate  the inner extremal coordinate (minimum z value)
     *  @param  outerCoordinate  the outer extremal coordinate (maximum z value)
     */
    void GetExtremalCoordinates(const pandora::CartesianPointList &pointList, pandora::CartesianVector &innerCoordinate, 
        pandora::CartesianVector &outerCoordinate) const;

    /**
     *  @brief Fit the start/end region of a given set of space points
     *
     *  @param  pointList  the input list of space points
     *  @param  vertexPosition  the input seed vertex position
     *  @param  vertexDirection  the input seed vertex direction 
     *  @param  fittedDirection  the output fitted vertex direction
     */
    void FitPoints(const pandora::CartesianPointList &pointList, const pandora::CartesianVector &vertexPosition,
        const pandora::CartesianVector &vertexDirection, pandora::CartesianVector &fittedDirection) const;


    bool                           m_isInitialized;    ///<
    pandora::CartesianVector       m_innerPosition;    ///<
    pandora::CartesianVector       m_innerDirection;   ///<
    pandora::CartesianVector       m_outerPosition;    ///<
    pandora::CartesianVector       m_outerDirection;   ///<
};

} // namespace lar_pandora

#endif // #ifndef PFPARTICLE_SEED_H
