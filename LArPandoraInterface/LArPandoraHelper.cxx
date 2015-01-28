/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.cxx
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */

#include "LArPandoraHelper.h"
#include "PFParticleSeed.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"
#include "ClusterFinder/ClusterCreator.h"

#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <limits>
#include <algorithm> // std::transform()
#include <iterator> // std::back_inserter()
#include <iostream>

namespace lar_pandora {

recob::Cluster LArPandoraHelper::BuildCluster(
  const int id, const std::vector<art::Ptr<recob::Hit>> &hitVector,
  cluster::ClusterParamsAlgBase& algo
)
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
    
    // Loop over vector of hits and calculate properties
    for (std::vector<art::Ptr<recob::Hit>>::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;

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
    
    // feed the algorithm with all the cluster hits;
    // usually this work is done by ClusterParamsImportWrapper<>,
    // but that is a static wrapper, and to use it we should make
    // this member function a template of the wrapped algorithm
    // (that would not be too bad anyway: users would barely notice,
    // but all this code should be moved into a header or a "tcc" file).
    // Instead, we are doing the hard work here:
    std::vector<recob::Hit const*> hits;
    std::transform(hitVector.begin(), hitVector.end(), std::back_inserter(hits),
      [](art::Ptr<recob::Hit> const& ptr) { return &*ptr; });
    algo.SetHits(hits);
    
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

recob::Track LArPandoraHelper::BuildTrack(const int id, const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints)
{
    mf::LogDebug("LArPandora") << "   Building Track [" << id << "], Number of Space Points = " << spacepoints.size() << std::endl;

    if (spacepoints.empty())
        throw cet::exception("LArPandora") << " LArPandoraHelper::BuildTrack --- No input hits were provided ";

    // Fill list of track properties
    std::vector<TVector3>               xyz;
    std::vector<TVector3>               pxpypz;
    std::vector< std::vector <double> > dQdx = std::vector< std::vector<double> >(0);
    std::vector<double>                 fitMomentum = std::vector<double>(2, util::kBogusD);

    // Loop over vector of space points
    PFParticleTrajectoryPointList trajectorypoints;
    PFParticleFitter::BuildTrajectoryPointList(spacepoints, trajectorypoints);

    if (trajectorypoints.empty())
        throw cet::exception("LArPandora") << " LArPandoraHelper::BuildTrack --- No trajectory points were found ";

    for (PFParticleTrajectoryPointList::const_iterator iter = trajectorypoints.begin(), iterEnd = trajectorypoints.end(); 
        iter != iterEnd; ++iter)
    {
        const PFParticleTrajectoryPoint &nextPoint = *iter;
        xyz.push_back(TVector3(nextPoint.m_position.GetX(), nextPoint.m_position.GetY(), nextPoint.m_position.GetZ()));
        pxpypz.push_back(TVector3(nextPoint.m_direction.GetX(), nextPoint.m_direction.GetY(), nextPoint.m_direction.GetZ())); 
    }

    // Return a new recob::Track object (of the Bezier flavour)
    return recob::Track(xyz, pxpypz, dQdx, fitMomentum, id);
}

} // namespace lar_pandora
