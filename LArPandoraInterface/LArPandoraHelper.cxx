/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.cxx
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */

// Framework includes
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"
#include "ClusterFinder/ClusterCreator.h"

// Pandora includes
#include "Objects/ParticleFlowObject.h"
#include "Objects/TrackState.h"
#include "Objects/Vertex.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

// Local LArrPandora includes
#include "LArPandoraInterface/LArPandoraHelper.h"
#include "LArPandoraInterface/PFParticleSeed.h"

// System includes
#include <limits>
#include <algorithm> // std::transform()
#include <iterator> // std::back_inserter()
#include <iostream>

namespace lar_pandora {

recob::Cluster LArPandoraHelper::BuildCluster(const int id, const std::vector<art::Ptr<recob::Hit>> &hitVector,
    const std::set<art::Ptr<recob::Hit>> &hitList, cluster::ClusterParamsAlgBase& algo) 
{
    mf::LogDebug("LArPandora") << "   Building Cluster [" << id << "], Number of hits = " << hitVector.size() << std::endl;

    if (hitVector.empty())
        throw cet::exception("LArPandora") << " LArPandoraHelper::BuildCluster --- No input hits were provided ";

    // Fill list of cluster properties
    geo::View_t view(geo::kUnknown);
    geo::PlaneID planeID;

    double startWire(+std::numeric_limits<float>::max()), sigmaStartWire(0.0);
    double startTime(+std::numeric_limits<float>::max()), sigmaStartTime(0.0);
    double endWire(-std::numeric_limits<float>::max()), sigmaEndWire(0.0);
    double endTime(-std::numeric_limits<float>::max()), sigmaEndTime(0.0);
    
    std::vector<recob::Hit const*> hits_for_params;
    hits_for_params.reserve(hitVector.size());
    
    for (std::vector<art::Ptr<recob::Hit>>::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        art::Ptr<recob::Hit> const& hit = *iter;
        
        const double thisWire(hit->WireID().Wire);
        const double thisWireSigma(0.5);
        const double thisTime(hit->PeakTime());
        const double thisTimeSigma(double(2.*hit->RMS()));
        const geo::View_t thisView(hit->View());
        const geo::PlaneID thisPlaneID(hit->WireID().planeID());

        if (geo::kUnknown == view)
        {
            view = thisView;
            planeID = thisPlaneID;
        }

        if (!(thisView == view && thisPlaneID == planeID))
        {
            throw cet::exception("LArPandora") << " LArPandoraHelper::BuildCluster --- Input hits have inconsistent plane IDs ";
        }
        
        hits_for_params.push_back(&*hit);
        
        if (hitList.count(hit))
            continue;

        if (thisWire < startWire || (thisWire == startWire && thisTime < startTime))
        {
            startWire = thisWire;
            sigmaStartWire = thisWireSigma;
            startTime = thisTime;
            sigmaStartTime = thisTimeSigma;
        }

        if (thisWire > endWire || (thisWire == endWire && thisTime > endTime))
        {
            endWire = thisWire;
            sigmaEndWire = thisWireSigma;
            endTime = thisTime;
            sigmaEndTime = thisTimeSigma;
        }

    }
    
    // feed the algorithm with all the cluster hits
    algo.SetHits(hits_for_params);
    
    // create the recob::Cluster directly in the vector
    return cluster::ClusterCreator(
      algo,                  // algo
      startWire,             // start_wire
      sigmaStartWire,        // sigma_start_wire
      startTime,             // start_tick
      sigmaStartTime,        // sigma_start_tick
      endWire,               // end_wire
      sigmaEndWire,          // sigma_end_wire
      endTime,               // end_tick
      sigmaEndTime,          // sigma_end_tick
      id,                    // ID
      view,                  // view
      planeID,               // plane
      recob::Cluster::Sentry // sentry
      ).move();
    
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
recob::Track LArPandoraHelper::BuildTrack(const int id, const pandora::ParticleFlowObject *const pPfo)
{
    mf::LogDebug("LArPandora") << "   Building Track [" << id << "], PdgCode = " << pPfo->GetParticleId() << std::endl;

    // Use sliding fits to calculate 3D trajectory points
    art::ServiceHandle<geo::Geometry> theGeometry;
    const float layerPitch((theGeometry->WirePitch(geo::kU) + theGeometry->WirePitch(geo::kV) + theGeometry->WirePitch(geo::kW))/3.f);
    const unsigned int layerWindow(20);

    std::vector<pandora::TrackState> trackStateVector;

    try
    {
        lar_content::LArPfoHelper::GetSlidingFitTrajectory(pPfo, layerWindow, layerPitch, trackStateVector);
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
    }

    return BuildTrack(id, trackStateVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Track LArPandoraHelper::BuildTrack(const int id, const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints, const bool isCosmic)
{
    mf::LogDebug("LArPandora") << "   Building Track [" << id << "], Number of Space Points = " << spacepoints.size() << std::endl;

    // Use linear regression to calculate 3D trajectory points
    std::vector<pandora::TrackState> trackStateVector;

    try
    {
        PFParticleFitter::GetLinearTrajectory(spacepoints, trackStateVector, isCosmic);
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
    }

    return BuildTrack(id, trackStateVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Track LArPandoraHelper::BuildTrack(const int id, const std::vector<pandora::TrackState> &trackStateVector)
{
    mf::LogDebug("LArPandora") << "   Building Track [" << id << "], Number of trajectory points = " << trackStateVector.size() << std::endl;

    if (trackStateVector.empty())
        throw cet::exception("LArPandora") << " LArPandoraHelper::BuildTrack --- No input trajectory points were provided ";

    // Fill list of track properties
    std::vector<TVector3>               xyz;
    std::vector<TVector3>               pxpypz;
    std::vector< std::vector <double> > dQdx = std::vector< std::vector<double> >(0);
    std::vector<double>                 fitMomentum = std::vector<double>(2, util::kBogusD);

    // Loop over trajectory points
    for (std::vector<pandora::TrackState>::const_iterator tIter = trackStateVector.begin(), tIterEnd = trackStateVector.end();
        tIter != tIterEnd; ++tIter)
    {
        const pandora::TrackState &nextPoint = *tIter;
        const pandora::CartesianVector position(nextPoint.GetPosition());
        const pandora::CartesianVector direction(nextPoint.GetMomentum().GetUnitVector());
        xyz.push_back(TVector3(position.GetX(), position.GetY(), position.GetZ()));
        pxpypz.push_back(TVector3(direction.GetX(), direction.GetY(), direction.GetZ()));
    }

    // Return a new recob::Track object (of the Bezier variety)
    return recob::Track(xyz, pxpypz, dQdx, fitMomentum, id);
}

} // namespace lar_pandora
