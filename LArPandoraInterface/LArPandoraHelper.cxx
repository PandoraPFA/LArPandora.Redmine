/**
 *  @file  larpandora/LArPandoraInterface/LArPandoraHelper.cxx
 *
 *  @brief helper function for LArPandoraInterface producer module
 *
 */

#include "LArPandoraHelper.h"
#include "PFParticleSeed.h"

#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

recob::Cluster LArPandoraHelper::BuildCluster(const int id, const std::vector<art::Ptr<recob::Hit>> &hitVector)
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
    double dQdW(0.0) /* , sigmadQdW(0.0) */;
    double dTdW(0.0) /* , sigmadTdW(0.0) */;

    double Sq(0.0), Sqx(0.0), Sqy(0.0), Sqxy(0.0), Sqxx(0.0);

    // Loop over vector of hits and calculate properties
    for (std::vector<art::Ptr<recob::Hit>>::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
    {
        const art::Ptr<recob::Hit> hit = *iter;

        const double thisCharge(hit->Integral());
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

        Sq   += thisCharge;
        Sqx  += thisCharge * thisWire;
        Sqy  += thisCharge * thisTime;
        Sqxx += thisCharge * thisWire * thisWire;
        Sqxy += thisCharge * thisWire * thisTime;
    }
    
    if (endWire >= startWire)
    {
        dQdW = Sq / (1.0 + endWire - startWire);
    //    sigmadQdW = 0.0;
    }
    else
    {
        throw cet::exception("LArPandora") << " LArPandora::BuildCluster --- Failed to find start and end wires ";
    }
    
    const double numerator(Sq * Sqxy - Sqx * Sqy);
    const double denominator(Sq * Sqxx - Sqx * Sqx);

    if (denominator > 0.0)
    {
        dTdW = numerator / denominator;
    //    sigmadTdW = 0.0;
    }

    // Return a new recob::Cluster object
    // FIXME fiiixmeeeeeeeee
    return recob::Cluster(
      startWire,             // start_wire
      sigmaStartWire,        // sigma_start_wire
      startTime,             // start_tick
      sigmaStartTime,        // sigma_start_tick
      dQdW,            // start_charge charge on the start wire FIXME (this should not be dQ/dW)
      std::atan(dTdW), // start_angle angle of the start of the cluster, in [-pi,pi] FIXME
      0.,             // start_opening opening angle at the start of the cluster  FIXME
      endWire,               // end_wire
      sigmaEndWire,          // sigma_end_wire
      endTime,               // end_tick
      sigmaEndTime,          // sigma_end_tick
      0.,             // end_charge charge on the end wire  FIXME
      0.,             // end_angle angle of the end of the cluster, in [-pi,pi]  FIXME
      0.,             // end_opening opening angle at the end of the cluster  FIXME
      0.,             // integral total charge from fitted shape of hits  FIXME
      0.,             // integral_stddev standard deviation of hit charge from fitted shape  FIXME
      0.,             // summedADC total charge from signal ADC of hits  FIXME
      0.,             // summedADC_stddev standard deviation of signal ADC of hits  FIXME
      hitVector.size(),      // n_hits
      0.,             // wires_over_hits wires covered by cluster, divided by number of hits  FIXME
      0.,             // width a measure of the cluster width  FIXME
      id,                    // ID
      view,                  // view
      planeID,               // plane
      recob::Cluster::Sentry // check that all the parameters are in
      );
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
