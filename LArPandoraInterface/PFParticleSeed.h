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
#include "Objects/CartesianVector.h"
#include "Objects/TrackState.h"

// ROOT includes
#include "TVector3.h"

namespace lar_pandora
{

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
     *  @param  axisDirection  the input seed vertex direction 
     *  @param  fittedDirection  the output fitted vertex direction
     */
    static void FitPoints(const pandora::CartesianPointList &pointList, const pandora::CartesianVector &axisDirection, 
        pandora::CartesianVector &fittedDirection);   

    /**
     *  @brief Build a vector of trajectory points from a vector of space points
     *
     *  @param  spacepoints  the input list of space points
     *  @param  trajectoryPoints  the output list of trajectory points
     *  @param  isCosmic  choice of vertex (true for highest vertex in Y, false for lowest vertex in Z) 
     */
    static void GetLinearTrajectory(const SpacePointVector &spacepoints, std::vector<pandora::TrackState> &trajectoryPoints,
        const bool isCosmic);
};

} // namespace lar_pandora

#endif // #ifndef PFPARTICLE_SEED_H
