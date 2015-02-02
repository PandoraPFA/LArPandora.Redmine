#include "PFParticleSeed.h"

#include <limits>

namespace lar_pandora {

PFParticleTrajectoryPoint::PFParticleTrajectoryPoint(const pandora::CartesianVector &position, const pandora::CartesianVector &direction,
    const float displacement) :
    m_position(position), m_direction(direction), m_displacement(displacement)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleSeed::PFParticleSeed(const SpacePointVector &spacepoints) : 
    m_isInitialized(false), 
    m_innerPosition(0.f, 0.f, 0.f),
    m_innerDirection(0.f, 0.f, 0.f),
    m_outerPosition(0.f, 0.f, 0.f),
    m_outerDirection(0.f, 0.f, 0.f)
{
    pandora::CartesianPointList pointList;
    PFParticleFitter::BuildPointList(spacepoints, pointList);
    this->Initialize(pointList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleSeed::PFParticleSeed(const pandora::CartesianPointList &pointList) : 
    m_isInitialized(false), 
    m_innerPosition(0.f, 0.f, 0.f),
    m_innerDirection(0.f, 0.f, 0.f),
    m_outerPosition(0.f, 0.f, 0.f),
    m_outerDirection(0.f, 0.f, 0.f)
{
    this->Initialize(pointList);
}
   
//------------------------------------------------------------------------------------------------------------------------------------------ 

void PFParticleSeed::Initialize(const pandora::CartesianPointList &pointList)
{
    try
    {
        PFParticleFitter::GetExtremalCoordinates(pointList, m_innerPosition, m_outerPosition);

        pandora::CartesianVector axisDirection((m_outerPosition - m_innerPosition).GetUnitVector());
        pandora::CartesianVector innerDirection(0.f, 0.f, 0.f);
        pandora::CartesianVector outerDirection(0.f, 0.f, 0.f);

        PFParticleFitter::FitPoints(pointList, m_innerPosition, axisDirection, innerDirection);
        PFParticleFitter::FitPoints(pointList, m_outerPosition, axisDirection, outerDirection);

        m_innerDirection = (innerDirection.GetDotProduct(axisDirection) > 0.f ? innerDirection : innerDirection * -1.f);
        m_outerDirection = (outerDirection.GetDotProduct(axisDirection) < 0.f ? outerDirection : outerDirection * -1.f);
        m_isInitialized = true;
    }
    catch (pandora::StatusCodeException& )
    {  
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

PFParticleSeed::~PFParticleSeed()
{

}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool PFParticleSeed::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
pandora::CartesianVector PFParticleSeed::GetInnerPosition() const 
{ 
    if (m_isInitialized)
        return m_innerPosition;  

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PFParticleSeed::GetInnerDirection() const 
{ 
    if (m_isInitialized)
        return m_innerDirection; 

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PFParticleSeed::GetOuterPosition() const 
{ 
    if (m_isInitialized)
        return m_outerPosition;  

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PFParticleSeed::GetOuterDirection() const 
{ 
    if (m_isInitialized)
        return m_outerDirection; 

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleFitter::BuildPointList(const SpacePointVector &spacepoints, pandora::CartesianPointList &pointList)
{
    for (SpacePointVector::const_iterator iter = spacepoints.begin(), iterEnd = spacepoints.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::SpacePoint> point = *iter; 
        pointList.push_back(pandora::CartesianVector(point->XYZ()[0], point->XYZ()[1], point->XYZ()[2]));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PFParticleFitter::BuildTrajectoryPointList(const SpacePointVector &spacepoints, PFParticleTrajectoryPointList &trajectoryPointList,
    const bool isCosmic)
{
    pandora::CartesianPointList pointList;
    PFParticleFitter::BuildPointList(spacepoints, pointList);

    if (pointList.size() > 1) // need at least two points for a trajectory
    {
        pandora::CartesianVector innerPosition(0.f, 0.f, 0.f);
        pandora::CartesianVector outerPosition(0.f, 0.f, 0.f);
        PFParticleFitter::GetExtremalCoordinates(pointList, innerPosition, outerPosition);

        const bool switchEnds(isCosmic && (outerPosition.GetY() > innerPosition.GetY()));
        const pandora::CartesianVector vtxPosition(switchEnds ? outerPosition : innerPosition);
        const pandora::CartesianVector endPosition(switchEnds ? innerPosition : outerPosition);
        const pandora::CartesianVector vtxDirection((endPosition - vtxPosition).GetUnitVector());

        // TODO: Fit the points here

        for (pandora::CartesianPointList::const_iterator iter = pointList.begin(), iterEnd = pointList.end(); iter != iterEnd; ++iter)
        {
            const pandora::CartesianVector &thisPosition = *iter;
            const float thisDisplacement(vtxDirection.GetDotProduct(thisPosition - vtxPosition));
            trajectoryPointList.push_back(PFParticleTrajectoryPoint(thisPosition, vtxDirection, thisDisplacement));
        }
    }

    std::sort(trajectoryPointList.begin(), trajectoryPointList.end());
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

void PFParticleFitter::FitPoints(const pandora::CartesianPointList &pointList, const pandora::CartesianVector &vertexPosition,
    const pandora::CartesianVector &vertexDirection, pandora::CartesianVector &fittedDirection)
{
    static const float m_cellSize(0.5f);
    static const float m_maxDisplacement(10.f);

    pandora::ClusterFitPointList fitPointList;

    for (pandora::CartesianPointList::const_iterator iter = pointList.begin(), iterEnd = pointList.end(); iter != iterEnd; ++iter)
    {
        const pandora::CartesianVector &thisPosition = *iter;
        const float thisDisplacement(vertexDirection.GetDotProduct(thisPosition - vertexPosition));

        if (thisDisplacement > m_maxDisplacement)
            continue;

        fitPointList.push_back(pandora::ClusterFitPoint(thisPosition, vertexDirection, m_cellSize, 1.f, 0));
    }
         
    pandora::ClusterFitResult fitResult;   
    pandora::ClusterFitHelper::FitPoints(fitPointList, fitResult);

    if (!fitResult.IsFitSuccessful())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    fittedDirection = fitResult.GetDirection();
}

} // namespace lar_pandora
