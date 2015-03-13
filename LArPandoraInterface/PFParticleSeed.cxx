#include "LArPandoraInterface/PFParticleSeed.h"

#include <limits>

namespace lar_pandora {

void PFParticleFitter::BuildPointList(const SpacePointVector &spacepoints, pandora::CartesianPointList &pointList)
{
    for (SpacePointVector::const_iterator iter = spacepoints.begin(), iterEnd = spacepoints.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::SpacePoint> point = *iter; 
        pointList.push_back(pandora::CartesianVector(point->XYZ()[0], point->XYZ()[1], point->XYZ()[2]));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleFitter::GetExtremalCoordinates(const pandora::CartesianPointList &pointList, pandora::CartesianVector &innerCoordinate, 
    pandora::CartesianVector &outerCoordinate)
{
    // Find extremal coordinates along X, Y and Z axes
    pandora::CartesianVector minPointX(+std::numeric_limits<float>::max(), 0.f, 0.f);
    pandora::CartesianVector maxPointX(-std::numeric_limits<float>::max(), 0.f, 0.f);

    pandora::CartesianVector minPointY(0.f, +std::numeric_limits<float>::max(), 0.f);
    pandora::CartesianVector maxPointY(0.f, -std::numeric_limits<float>::max(), 0.f);

    pandora::CartesianVector minPointZ(0.f, 0.f, +std::numeric_limits<float>::max());
    pandora::CartesianVector maxPointZ(0.f, 0.f, -std::numeric_limits<float>::max());

    for (pandora::CartesianPointList::const_iterator iter = pointList.begin(), iterEnd = pointList.end(); iter != iterEnd; ++iter)
    {
        const pandora::CartesianVector &thisPoint = *iter;

        if (thisPoint.GetX() < minPointX.GetX()) 
            minPointX = thisPoint;

        if (thisPoint.GetX() > maxPointX.GetX()) 
            maxPointX = thisPoint;

        if (thisPoint.GetY() < minPointY.GetY()) 
            minPointY = thisPoint;

        if (thisPoint.GetY() > maxPointY.GetY()) 
            maxPointY = thisPoint;

        if (thisPoint.GetZ() < minPointZ.GetZ()) 
            minPointZ = thisPoint;

        if (thisPoint.GetZ() > maxPointZ.GetZ()) 
            maxPointZ = thisPoint;
    }

    if ((minPointX.GetX() > maxPointX.GetX()) || (minPointY.GetY() > maxPointY.GetY()) || (minPointZ.GetZ() > maxPointZ.GetZ()))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    // Find the pair of coordinates with the maximum separation
    pandora::CartesianPointList candidateList;
    candidateList.push_back(minPointX);
    candidateList.push_back(maxPointX);
    candidateList.push_back(minPointY);
    candidateList.push_back(maxPointY);
    candidateList.push_back(minPointZ);
    candidateList.push_back(maxPointZ);

    bool foundExtremalCoordinates(false);
    float maxDistanceSquared(std::numeric_limits<float>::epsilon());

    for (pandora::CartesianPointList::const_iterator iter1 = candidateList.begin(), iterEnd1 = candidateList.end(); iter1 != iterEnd1; ++iter1)
    {
        const pandora::CartesianVector &point1 = *iter1;

        for (pandora::CartesianPointList::const_iterator iter2 = iter1, iterEnd2 = candidateList.end(); iter2 != iterEnd2; ++iter2)
        {
            const pandora::CartesianVector &point2 = *iter2;

            const float thisDistanceSquared((point1 - point2).GetMagnitudeSquared());

            if (thisDistanceSquared > maxDistanceSquared)
            {
                maxDistanceSquared = thisDistanceSquared;
                innerCoordinate = ((point1.GetZ() < point2.GetZ()) ? point1 : point2);
                outerCoordinate = ((point1.GetZ() < point2.GetZ()) ? point2 : point1);
                foundExtremalCoordinates = true;
            }
        }
    }

    if (!foundExtremalCoordinates)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleFitter::FitPoints(const pandora::CartesianPointList &pointList, const pandora::CartesianVector &axisDirection, 
    pandora::CartesianVector &fittedDirection)
{
    static const float m_cellSize(0.5f);

    pandora::ClusterFitPointList fitPointList;

    for (pandora::CartesianPointList::const_iterator iter = pointList.begin(), iterEnd = pointList.end(); iter != iterEnd; ++iter)
    {
        const pandora::CartesianVector &thisPosition = *iter;
        fitPointList.push_back(pandora::ClusterFitPoint(thisPosition, axisDirection, m_cellSize, 1.f, 0));
    }
         
    pandora::ClusterFitResult fitResult;   
    pandora::ClusterFitHelper::FitPoints(fitPointList, fitResult);

    if (!fitResult.IsFitSuccessful())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    fittedDirection = fitResult.GetDirection();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleFitter::GetLinearTrajectory(const SpacePointVector &spacepoints, std::vector<pandora::TrackState> &trajectoryPoints,
    const bool isCosmic)
{
    try
    {
        // Need at least two points for a linear fit 
        if (spacepoints.size() < 2)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

        // Convert space points into cartesian vector points
        pandora::CartesianPointList pointList;
        PFParticleFitter::BuildPointList(spacepoints, pointList);

        // Find inner and outer positions
        pandora::CartesianVector innerPosition(0.f, 0.f, 0.f);
        pandora::CartesianVector outerPosition(0.f, 0.f, 0.f);
        PFParticleFitter::GetExtremalCoordinates(pointList, innerPosition, outerPosition);
        
        // Choose vertex direction
        const bool switchEnds(isCosmic && (outerPosition.GetY() > innerPosition.GetY()));
        const pandora::CartesianVector vtxPosition(switchEnds ? outerPosition : innerPosition);

        // Perform linear regression and choose vertex direction
        const pandora::CartesianVector seedDirection((outerPosition - innerPosition).GetUnitVector());
        pandora::CartesianVector axisDirection(0.f, 0.f, 0.f);        
        PFParticleFitter::FitPoints(pointList, seedDirection, axisDirection);

        const float firstCorrection((seedDirection.GetDotProduct(axisDirection) < 0.f) ? -1.f : +1.f);
        const float secondCorrection(switchEnds ? -1.f : +1.f);
        const pandora::CartesianVector vtxDirection(axisDirection * firstCorrection * secondCorrection);

        // Create ordered trajectory list
        typedef std::map< const float, const pandora::TrackState > ThreeDTrajectoryMap;
        ThreeDTrajectoryMap trajectoryMap;

        for (pandora::CartesianPointList::const_iterator iter = pointList.begin(), iterEnd = pointList.end(); iter != iterEnd; ++iter)
        {
            const pandora::CartesianVector &thisPosition = *iter;
            const float displacement(vtxDirection.GetDotProduct(thisPosition - vtxPosition));    
            trajectoryMap.insert(std::pair<const float, const pandora::TrackState>(displacement, pandora::TrackState(thisPosition, vtxDirection)));
	}
     
        for (ThreeDTrajectoryMap::const_iterator tIter = trajectoryMap.begin(), tIterEnd = trajectoryMap.end(); tIter != tIterEnd; ++tIter)
	{
            trajectoryPoints.push_back(tIter->second);
	}
    }
    catch (pandora::StatusCodeException& )
    {  
    }
}

} // namespace lar_pandora
