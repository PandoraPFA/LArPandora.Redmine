///////////////////////////////////////////////////////////////////////
///
/// \file   SpacePointService.cxx
///
/// \brief  Service for generating space points from hits.
///
/// \author H. Greenlee 
///
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "Simulation/SimChannel.h"
#include "art/Framework/Principal/View.h"
#include "Utilities/DetectorProperties.h"
#include "MCCheater/BackTracker.h"
#include "TH1F.h"

namespace {

  // Function classes for sorting sim::IDEs according to track id.

  //class IDELess {
  //public:
  //  bool operator()(const sim::IDE& p1, const sim::IDE& p2) {
  //    bool result = p1.trackID < p2.trackID;
  //    return result;
  //  }
  //};

  //class IDEEqual {
  //public:
  //  bool operator()(const sim::IDE& p1, const sim::IDE& p2) {
  //    bool result = p1.trackID == p2.trackID;
  //    return result;
  //  }
  //};

}

//----------------------------------------------------------------------
// Constructor.
//
trkf::SpacePointService::SpacePointService(const fhicl::ParameterSet& pset,
					   art::ActivityRegistry& reg) :
  fUseMC(false),
  fMaxDT(0.),
  fMaxS(0.),
  fTimeOffsetU(0.),
  fTimeOffsetV(0.),
  fTimeOffsetW(0.),
  fMinViews(1000),
  fEnableU(false),
  fEnableV(false),
  fEnableW(false),
  fFilter(false)
{
  reconfigure(pset);
}

//----------------------------------------------------------------------
// Destructor.
//
trkf::SpacePointService::~SpacePointService()
{
}

//----------------------------------------------------------------------
// Update configuration parameters.
//
void trkf::SpacePointService::reconfigure(const fhicl::ParameterSet& pset)
{
  // Get configuration parameters.

  fUseMC = pset.get<bool>("UseMC", false);
  fMaxDT = pset.get<double>("MaxDT", 0.);
  fMaxS = pset.get<double>("MaxS", 0.);

  fTimeOffsetU = pset.get<double>("TimeOffsetU", 0.);
  fTimeOffsetV = pset.get<double>("TimeOffsetV", 0.);
  fTimeOffsetW = pset.get<double>("TimeOffsetW", 0.);

  fMinViews = pset.get<int>("MinViews", 1000);

  fEnableU = pset.get<bool>("EnableU", false);
  fEnableV = pset.get<bool>("EnableV", false);
  fEnableW = pset.get<bool>("EnableW", false);
  fFilter = pset.get<bool>("Filter", false);

  // Report.

  mf::LogInfo("SpacePointService") 
    << "SpacePointService configured with the following parameters:\n"
    << "  UseMC = " << fUseMC << "\n"
    << "  MaxDT = " << fMaxDT << "\n"
    << "  MaxS = " << fMaxS << "\n"
    << "  TimeOffsetU = " << fTimeOffsetU << "\n"
    << "  TimeOffsetV = " << fTimeOffsetV << "\n"
    << "  TimeOffsetW = " << fTimeOffsetW << "\n" 
    << "  MinViews = " << fMinViews << "\n"
    << "  EnableU = " << fEnableU << "\n"
    << "  EnableV = " << fEnableV << "\n"
    << "  EnableW = " << fEnableW << "\n"
    << "  Filter = " << fFilter;
}

//----------------------------------------------------------------------
// Print geometry and properties constants.
//
void trkf::SpacePointService::update() const
{
  // Generate info report on first call only.

  static bool first = true;
  bool report = first;
  first = false;

  // Get services.

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> larprop;

  // Calculate and print geometry information.

  mf::LogInfo log("SpacePointService");
  if(report)
    log << "Updating geometry constants.\n";

  // Update detector properties.

  double samplingRate = detprop->SamplingRate();
  double triggerOffset = detprop->TriggerOffset();
  if(report) {
    log << "\nDetector properties:\n"
	<< "  Sampling Rate = " << samplingRate << " ns/tick\n"
	<< "  Trigger offset = " << triggerOffset << " ticks\n";
  }

  // Update LArProperties.

  double efield = larprop->Efield();
  double temperature = larprop->Temperature();
  double driftVelocity = larprop->DriftVelocity(efield, temperature);
  double timePitch = 0.001 * driftVelocity * samplingRate;
  if(report) {
    log << "\nLAr propertoes:\n"
	<< "  E field = " << efield << " kV/cm\n"
	<< "  Temperature = " << temperature << " K\n"
	<< "  Drift velocity = " << driftVelocity << " cm/us\n"
	<< "  Time pitch = " << timePitch << " cm/tick\n";
  }

  // Get time offsets.

  std::vector<std::vector<double> > timeOffset;
  fillTimeOffset(timeOffset);

  // Loop over TPCs.

  int ntpc = geom->NTPC();

  for(int tpc = 0; tpc < ntpc; ++tpc) {
    const geo::TPCGeo& tpcgeom = geom->TPC(tpc);

    // Loop over planes.

    int nplane = tpcgeom.Nplanes();

    for(int plane = 0; plane < nplane; ++plane) {
      const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);

      // Fill view-dependent quantities.

      geo::View_t view = pgeom.View();
      std::string viewname = "?";
      if(view == geo::kU) {
	viewname = "U";
      }
      else if(view == geo::kV) {
	viewname = "V";
      }
      else if(view == geo::kW) {
	viewname = "W";
      }
      else
	throw cet::exception("SpacePointService") << "Bad view = " 
						  << view << "\n";

      std::string sigtypename = "?";
      geo::SigType_t sigtype = pgeom.SignalType();
      if(sigtype == geo::kInduction)
	sigtypename = "Induction";
      else if(sigtype == geo::kCollection)
	sigtypename = "Collection";
      else
	throw cet::exception("SpacePointService") << "Bad signal type = " 
						  << sigtype << "\n";

      std::string orientname = "?";
      geo::Orient_t orient = pgeom.Orientation();
      if(orient == geo::kVertical)
	orientname = "Vertical";
      else if(orient == geo::kHorizontal)
	orientname = "Horizontal";
      else
	throw cet::exception("SpacePointService") << "Bad orientation = " 
						  << orient << "\n";

      if(report) {
	const double* xyz = tpcgeom.PlaneLocation(plane);
	log << "\nTPC, Plane: " << tpc << ", " << plane << "\n"
	    << "  View: " << viewname << "\n"
	    << "  SignalType: " << sigtypename << "\n"
	    << "  Orientation: " << orientname << "\n"
	    << "  Plane location: " << xyz[0] << "\n"
	    << "  Plane pitch: " << tpcgeom.Plane0Pitch(plane) << "\n"
	    << "  Wire angle: " << tpcgeom.Plane(plane).Wire(0).ThetaZ() << "\n"
	    << "  Wire pitch: " << tpcgeom.WirePitch() << "\n"
	    << "  Time offset: " << timeOffset[tpc][plane] << "\n";
      }

      if(orient != geo::kVertical)
	throw cet::exception("SpacePointService") 
	  << "Horizontal wire geometry not implemented.\n";
    }
  }
}

//----------------------------------------------------------------------
// Calculate time offsets.
// Results stored in nested vector indexed by [tpc][plane]
void trkf::SpacePointService::fillTimeOffset(std::vector<std::vector<double> >& timeOffset) const
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> larprop;

  // Clear result.

  timeOffset.clear();

  // Get properties needed to calculate time offsets.

  double samplingRate = detprop->SamplingRate();
  double triggerOffset = detprop->TriggerOffset();
  double efield = larprop->Efield();
  double temperature = larprop->Temperature();
  double driftVelocity = larprop->DriftVelocity(efield, temperature);
  double timePitch = 0.001 * driftVelocity * samplingRate;

  // Loop over TPCs.

  int ntpc = geom->NTPC();
  timeOffset.resize(ntpc);

  for(int tpc = 0; tpc < ntpc; ++tpc) {
    const geo::TPCGeo& tpcgeom = geom->TPC(tpc);

    // Loop over planes.

    int nplane = tpcgeom.Nplanes();
    timeOffset[tpc].resize(nplane, 0.);

    for(int plane = 0; plane < nplane; ++plane) {
      const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);

      // Calculate geometric time offset.

      const double* xyz = tpcgeom.PlaneLocation(0);
      timeOffset[tpc][plane] =
	(-xyz[0] + tpcgeom.Plane0Pitch(plane)) / timePitch + triggerOffset;

      // Add view-dependent time offset.

      geo::View_t view = pgeom.View();
      if(view == geo::kU)
	timeOffset[tpc][plane] += fTimeOffsetU;
      else if(view == geo::kV)
	timeOffset[tpc][plane] += fTimeOffsetV;
      else if(view == geo::kW)
	timeOffset[tpc][plane] += fTimeOffsetW;
      else
	throw cet::exception("SpacePointService") << "Bad view = " 
						  << view << "\n";
    }
  }
}



//----------------------------------------------------------------------
// Get corrected time for the specified hit.
double trkf::SpacePointService::correctedTime(const recob::Hit& hit,
					      const std::vector<std::vector<double> >& timeOffset) const
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geom;

  // Get tpc, plane.

  unsigned short channel = hit.Channel();
  unsigned int tpc, plane, wire;
  geom->ChannelToWire(channel, tpc, plane, wire);

  // Correct time for trigger offset and plane-dependent time offsets.

  double t = hit.PeakTime() - timeOffset[tpc][plane];

  return t;
}

//----------------------------------------------------------------------
// Spatial separation of hits (zero if two or fewer).
double trkf::SpacePointService::separation(const art::PtrVector<recob::Hit>& hits) const
{
  // Get geometry service.

  art::ServiceHandle<geo::Geometry> geom;

  // Trivial case - fewer than three hits.

  if(hits.size() < 3)
    return 0.;

  // Error case - more than three hits.

  if(hits.size() > 3) {
    mf::LogError("SpacePointService") << "Method separation called with more than three htis.";
      return 0.;
  }

  // Got exactly three hits.

  assert(hits.size() == 3);

  // Calculate angles and distance of each hit from origin.

  double dist[3] = {0., 0., 0.};
  double sinth[3] = {0., 0., 0.};
  double costh[3] = {0., 0., 0.};
  unsigned int tpcs[3];
  unsigned int planes[3];

  for(int i=0; i<3; ++i) {

    // Get tpc, plane, wire.

    const recob::Hit& hit = *(hits[i]);
    unsigned short channel = hit.Channel();
    unsigned int tpc, plane, wire;
    const geo::WireGeo& wgeom = geom->ChannelToWire(channel, tpc, plane, wire);
    tpcs[i] = tpc;
    planes[i] = plane;

    // Check tpc and plane errors.

    for(int j=0; j<i; ++j) {

      if(tpcs[j] != tpc) {
	mf::LogError("SpacePointService") << "Method separation called with hits from multiple tpcs..";
	return 0.;
      }

      if(planes[j] == plane) {
	mf::LogError("SpacePointService") << "Method separation called with hits from the same plane..";
	return 0.;
      }
    }

    // Get angles and distance of wire.

    double hl = wgeom.HalfL();
    double xyz[3];
    double xyz1[3];
    wgeom.GetCenter(xyz);
    wgeom.GetCenter(xyz1, hl);
    double s = (xyz1[1] - xyz[1]) / hl;
    double c = (xyz1[2] - xyz[2]) / hl;
    sinth[plane] = s;
    costh[plane] = c;
    dist[plane] = xyz[2] * s - xyz[1] * c;
  }

  double S = ((sinth[1] * costh[2] - costh[1] * sinth[2]) * dist[0] 
	      +(sinth[2] * costh[0] - costh[2] * sinth[0]) * dist[1] 
	      +(sinth[0] * costh[1] - costh[0] * sinth[1]) * dist[2]);
  return S;
}

//----------------------------------------------------------------------
// Check hits for compatibility.
// Check hits pairwise for different views and maximum time difference.
// Check three hits for spatial compatibility.
bool trkf::SpacePointService::compatible(const art::PtrVector<recob::Hit>& hits,
					 const std::vector<std::vector<double> >& timeOffset,
					 double maxDT, double maxS) const
{
  // Get geometry service.

  art::ServiceHandle<geo::Geometry> geom;

  // Get cuts.

  if(maxDT == 0.)
    maxDT = fMaxDT;
  if(maxS == 0.)
    maxS = fMaxS;  

  int nhits = hits.size();

  // Fewer than two or more than three hits can never be compatible.

  bool result = nhits >= 2 && nhits <= 3;
  bool mc_ok = true;
  unsigned int tpc = 0;

  if(result) {

    // First do pairwise tests.
    // Do double loop over hits.

    for(int ihit1 = 0; result && ihit1 < nhits-1; ++ihit1) {
      const recob::Hit& hit1 = *(hits[ihit1]);
      unsigned short channel1 = hit1.Channel();
      unsigned int tpc1, plane1, wire1;
      geom->ChannelToWire(channel1, tpc1, plane1, wire1);
      geo::View_t view1 = hit1.View();
      double t1 = hit1.PeakTime() - timeOffset[tpc1][plane1];

      // If using mc information, get a collection of track ids for hit 1.

      const HitMCInfo& mcinfo1 = fHitMCMap[&hit1];
      const std::vector<int>& tid1 = mcinfo1.trackIDs;
      bool only_neg1 = tid1.size() > 0 && tid1.back() < 0;

      // Loop over second hit.

      for(int ihit2 = ihit1+1; result && ihit2 < nhits; ++ihit2) {
	const recob::Hit& hit2 = *(hits[ihit2]);
	unsigned short channel2 = hit2.Channel();
	unsigned int tpc2, plane2, wire2;
	geom->ChannelToWire(channel2, tpc2, plane2, wire2);
	geo::View_t view2 = hit2.View();

	// Test for same tpc and different views.

	result = result && tpc1 == tpc2 && view1 != view2;
	if(result) {

	  // Remember which tpc we are in.

	  tpc = tpc1;

	  double t2 = hit2.PeakTime() - timeOffset[tpc2][plane2];
    
	  // Test maximum time difference.

	  result = result && std::abs(t1-t2) <= maxDT;

	  // Test mc truth.

	  if(result && fUseMC) {

	    // Test whether hits have a common parent track id.

	    const HitMCInfo& mcinfo2 = fHitMCMap[&hit2];
	    std::vector<int> tid2 = mcinfo2.trackIDs;
	    bool only_neg2 = tid2.size() > 0 && tid2.back() < 0;
	    std::vector<int>::iterator it =
	      std::set_intersection(tid1.begin(), tid1.end(),
				    tid2.begin(), tid2.end(),
				    tid2.begin());
	    tid2.resize(it - tid2.begin());

	    // Hits are compatible if they have parents in common.
	    // If the only parent id in common is negative (-999),
	    // then hits are compatible only if both hits have only
	    // negative parent tracks.

	    bool only_neg3 = tid2.size() > 0 && tid2.back() < 0;
	    mc_ok = tid2.size() > 0 && 
	      (!only_neg3 || (only_neg1 && only_neg2));
	    result = result && mc_ok;

	    // If we are still OK, check that either hit is
	    // the nearest neighbor of the other.

	    if(result) {
	      result = mcinfo1.pchit[plane2] == &hit2 || 
		mcinfo2.pchit[plane1] == &hit1;
	    }
	  }
	}
      }
    }

    // If there are exactly three hits, and they pass pairwise tests, check
    // for spatial compatibility.

    if(result && nhits == 3) {

      // Loop over hits.

      double dist[3] = {0., 0., 0.};
      double sinth[3] = {0., 0., 0.};
      double costh[3] = {0., 0., 0.};
      geo::View_t view[3];

      for(int i=0; i<3; ++i) {

	// Get tpc, plane, wire.

	const recob::Hit& hit = *(hits[i]);
	unsigned short channel = hit.Channel();
	unsigned int tpc0, plane, wire;
	const geo::WireGeo& wgeom = geom->ChannelToWire(channel, tpc0, plane, wire);
	assert(tpc0 == tpc);
	view[i] = hit.View();

	// Get angles and distance of wire.

	double hl = wgeom.HalfL();
	double xyz[3];
	double xyz1[3];
	wgeom.GetCenter(xyz);
	wgeom.GetCenter(xyz1, hl);
	double s  = (xyz1[1] - xyz[1]) / hl;
	double c = (xyz1[2] - xyz[2]) / hl;
	sinth[plane] = s;
	costh[plane] = c;
	dist[plane] = xyz[2] * s - xyz[1] * c;
      }

      // Do space cut.

      double S = ((sinth[1] * costh[2] - costh[1] * sinth[2]) * dist[0] 
		  +(sinth[2] * costh[0] - costh[2] * sinth[0]) * dist[1] 
		  +(sinth[0] * costh[1] - costh[0] * sinth[1]) * dist[2]);

      result = result && std::abs(S) < maxS;
    }
  }

  // Done.
    
  return result;
}

//----------------------------------------------------------------------
// Fill one space point using a colleciton of hits.
// Assume points have already been tested for compatibility.
// Returned value is a time goodness of fit (smaller is better, like chisquare).
//
double trkf::SpacePointService::fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
					       const std::vector<std::vector<double> >& timeOffset,
					       recob::SpacePoint& spt) const
{
  // Get services.

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArProperties> larprop;

  // Calculate time pitch.

  double efield = larprop->Efield();
  double temperature = larprop->Temperature();
  double driftVelocity = larprop->DriftVelocity(efield, temperature);
  double samplingRate = detprop->SamplingRate();
  double timePitch = 0.001 * driftVelocity * samplingRate;

  // Store hits in SpacePoint.

  spt = recob::SpacePoint(recob::SpacePoint(hits));
  int nhits = hits.size();

  // Calculate position.

  double xyz[3] = {0., 0., 0.};

  // Calculate x using drift times.
  // Loop over all hits and calculate the weighted average drift time.
  // Also calculate time variance, minimum and maximum time.

  double sumt2w = 0.;
  double sumtw = 0.;
  double sumw = 0.;

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
      ihit != hits.end(); ++ihit) {

    const recob::Hit& hit = **ihit;
    unsigned short channel = hit.Channel();
    unsigned int tpc, plane, wire;
    geom->ChannelToWire(channel, tpc, plane, wire);

    // Correct time for trigger offset and view-dependent time offsets.
    // Assume time error is proportional to (end time - start time).

    double t0 = timeOffset[tpc][plane];
    double t = hit.PeakTime() - t0;
    double t1 = hit.StartTime() - t0;
    double t2 = hit.EndTime() - t0;
    double et = t2 - t1;
    double w = 1./(et*et);

    sumt2w += w*t*t;
    sumtw += w*t;
    sumw += w;
  }

  double drift_time = 0.;
  double var = 0.;
  if(sumw != 0.) {
    drift_time = sumtw / sumw;
    var = sumt2w / sumw - drift_time * drift_time;
  }
  xyz[0] = drift_time * timePitch;
  var *= timePitch;

  // Calculate y, z using wires (need at least two hits).

  if(nhits >= 2) {

    // Calculate y and z by chisquare minimization of wire coordinates.

    double sus = 0.;   // sum u_i sin_th_i
    double suc = 0.;   // sum u_i cos_th_i
    double sc2 = 0.;   // sum cos2_th_i
    double ss2 = 0.;   // sum sin2_th_i
    double ssc = 0.;   // sum sin_th_i cos_th_i

    // Loop over points.

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;
      unsigned short channel = hit.Channel();
      unsigned int tpc, plane, wire;
      const geo::WireGeo& wgeom = geom->ChannelToWire(channel, tpc, plane, wire);

      // Calculate angle and wire coordinate in this view.
    
      double hl = wgeom.HalfL();
      double xyz[3];
      double xyz1[3];
      wgeom.GetCenter(xyz);
      wgeom.GetCenter(xyz1, hl);
      double s  = (xyz1[1] - xyz[1]) / hl;
      double c = (xyz1[2] - xyz[2]) / hl;
      double u = xyz[2] * s - xyz[1] * c;

      // Summations

      sus += u*s;
      suc += u*c;
      sc2 += c*c;
      ss2 += s*s;
      ssc += s*c;
    }

    // Calculate y,z

    double denom = sc2 * ss2 - ssc * ssc;
    if(denom != 0.) {
      xyz[1] = (-suc * ss2 + sus * ssc) / denom;
      xyz[2] = (sus * sc2 - suc * ssc) / denom;
    }

    // Set coordintates in space point.

    spt.SetXYZ(xyz);
  }
  return var;
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits (non-mc-truth version).
//
void trkf::SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					      std::vector<recob::SpacePoint>& spts,
					      double maxDT, double maxS) const
{
  std::vector<const sim::SimChannel*> empty;
  makeSpacePoints(hits, spts, empty, maxDT, maxS);
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits (mc truth version).
//
void trkf::SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					      std::vector<recob::SpacePoint>& spts,
					      const std::vector<const sim::SimChannel*>& simchans,
					      double maxDT, double maxS) const
{
  // Get cuts.

  if(maxDT == 0.)
    maxDT = fMaxDT;
  if(maxS == 0.)
    maxS = fMaxS;  

  // Get geometry service.

  art::ServiceHandle<geo::Geometry> geom;

  // Print diagnostic information.

  update();

  // Get time offsets.

  std::vector<std::vector<double> > timeOffset;
  fillTimeOffset(timeOffset);

  // First make result vector is empty.

  spts.erase(spts.begin(), spts.end());

  // Statistics.

  int n2 = 0;  // Number of two-hit space points.
  int n3 = 0;  // Number of three-hit space points.
  int n2filt = 0;  // Number of two-hit space points after filtering.
  int n3filt = 0;  // Number of three-hit space pointe after filtering.

  // If fUseMC is true, verify that channels are sorted by channel number.

  if(fUseMC) {

    unsigned int nsc = simchans.size();
    for(unsigned int isc = 0; isc < nsc; ++isc) {
      const sim::SimChannel* psc = simchans[isc];
      if(psc != 0 && isc != psc->Channel())
	throw cet::exception("SpacePointService") << "MC channels not sorted.\n";
    }
  }

  // Sort hits into maps indexed by [tpc][plane][wire].
  // If using mc information, also generate maps of sim::IDEs and mc 
  // position indexed by hit.

  std::vector<std::vector<std::map<unsigned int, art::Ptr<recob::Hit> > > > hitmap;
  fHitMCMap.clear();

  unsigned int ntpc = geom->NTPC();
  hitmap.resize(ntpc);
  for(unsigned int tpc = 0; tpc < ntpc; ++tpc) {
    int nplane = geom->Nplanes(tpc);
    hitmap[tpc].resize(nplane);
  }

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
      ihit != hits.end(); ++ihit) {
    const art::Ptr<recob::Hit>& phit = *ihit;
    geo::View_t view = phit->View();
    if((view == geo::kU && fEnableU) ||
       (view == geo::kV && fEnableV) ||
       (view == geo::kW && fEnableW)) {

      unsigned short channel = phit->Channel();
      unsigned int tpc, plane, wire;
      geom->ChannelToWire(channel, tpc, plane, wire);
      hitmap[tpc][plane][wire] = phit;
    }
  }

  // Fill mc information, including IDEs and closest neighbors
  // of each hit.

  if(fUseMC) {

    // First loop over hits and fill track ids and mc position.

    for(unsigned int tpc = 0; tpc < ntpc; ++tpc) {
      int nplane = geom->Nplanes(tpc);
      for(int plane = 0; plane < nplane; ++plane) {
	for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator ihit = hitmap[tpc][plane].begin();
	    ihit != hitmap[tpc][plane].end(); ++ihit) {
	  const art::Ptr<recob::Hit>& phit = ihit->second;
	  const recob::Hit& hit = *phit;
	  HitMCInfo& mcinfo = fHitMCMap[&hit];   // Default HitMCInfo.

	  // Fill default nearest neighbor information (i.e. none).

	  mcinfo.pchit.resize(nplane, 0);
	  mcinfo.dist2.resize(nplane, 1.e20);

	  // Get sim::IDEs for this hit.

	  std::vector<sim::IDE> ides;
	  cheat::BackTracker::HitToSimIDEs(*simchans[hit.Channel()], phit, ides);

	  // Get sorted track ids. for this hit.

	  mcinfo.trackIDs.reserve(ides.size());
	  for(std::vector<sim::IDE>::const_iterator i = ides.begin();
	      i != ides.end(); ++i)
	    mcinfo.trackIDs.push_back(i->trackID);
	  sort(mcinfo.trackIDs.begin(), mcinfo.trackIDs.end());

	  // Get position of ionization for this hit.

	  mcinfo.xyz = cheat::BackTracker::SimIDEsToXYZ(ides);
	}
      }
    }

    // Loop over hits again and fill nearest neighbor information for real.

    for(unsigned int tpc = 0; tpc < ntpc; ++tpc) {
      int nplane = geom->Nplanes(tpc);
      for(int plane = 0; plane < nplane; ++plane) {
	for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator ihit = hitmap[tpc][plane].begin();
	    ihit != hitmap[tpc][plane].end(); ++ihit) {
	  const art::Ptr<recob::Hit>& phit = ihit->second;
	  const recob::Hit& hit = *phit;
	  HitMCInfo& mcinfo = fHitMCMap[&hit];

	  // Fill nearest neighbor information for this hit.

	  for(int plane2 = 0; plane2 < nplane; ++plane2) {
	    for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator jhit = hitmap[tpc][plane2].begin();
		jhit != hitmap[tpc][plane2].end(); ++jhit) {
	      const art::Ptr<recob::Hit>& phit2 = jhit->second;
	      const recob::Hit& hit2 = *phit2;
	      const HitMCInfo& mcinfo2 = fHitMCMap[&hit2];

	      assert(mcinfo.xyz.size() == 3);
	      assert(mcinfo2.xyz.size() == 3);
	      double dx = mcinfo.xyz[0] - mcinfo2.xyz[0];
	      double dy = mcinfo.xyz[1] - mcinfo2.xyz[1];
	      double dz = mcinfo.xyz[2] - mcinfo2.xyz[2];
	      double dist2 = dx*dx + dy*dy + dz*dz;
	      if(dist2 < mcinfo.dist2[plane2]) {
		mcinfo.dist2[plane2] = dist2;
		mcinfo.pchit[plane2] = &hit2;
	      }
	    }
	  }
	}
      }
    }
  }

  mf::LogDebug debug("SpacePointService");
  debug << "Total hits = " << hits.size() << "\n\n";

  for(unsigned int tpc = 0; tpc < ntpc; ++tpc) {
    int nplane = hitmap[tpc].size();
    for(int plane = 0; plane < nplane; ++plane) {
      debug << "TPC, Plane: " << tpc << ", " << plane 
	    << ", hits = " << hitmap[tpc][plane].size() << "\n";
    }
  }

  // Loop over TPCs.

  for(unsigned int tpc = 0; tpc < ntpc; ++tpc) {

    // Make empty multimap from hit pointer on most-populated plane to space points 
    // that include that hit (used for filtering).

    typedef const recob::Hit* sptkey_type;
    std::multimap<sptkey_type, SpacePointX> sptmap;
    std::set<sptkey_type> sptkeys;              // Keys of multimap.

    // Sort maps in increasing order of number of hits.
    // This is so that we can do the outer loops over hits 
    // over the views with fewer hits.

    int nplane = hitmap[tpc].size();
    std::vector<int> index(nplane);

    for(int i=0; i<nplane; ++i)
      index[i] = i;

    for(int i=0; i<nplane-1; ++i) {
      for(int j=i+1; j<nplane; ++j) {
	if(hitmap[tpc][index[i]].size() > hitmap[tpc][index[j]].size()) {
	  int temp = index[i];
	  index[i] = index[j];
	  index[j] = temp;
	}
      }
    }

    // If two-view space points are allowed, make a double loop
    // over hits and produce space points for compatible hit-pairs.

    if(fMinViews <= 2) {

      // Loop over pairs of views.

      for(int i=0; i<nplane-1; ++i) {
	unsigned int plane1 = index[i];

	for(int j=i+1; j<nplane; ++j) {
	  unsigned int plane2 = index[j];

	  // Get angle, pitch, and offset of plane2 wires.

	  const geo::WireGeo& wgeo2 = geom->Plane(plane2, tpc).Wire(0);
	  double hl2 = wgeo2.HalfL();
	  double xyz21[3];
	  double xyz22[3];
	  wgeo2.GetCenter(xyz21, -hl2);
	  wgeo2.GetCenter(xyz22, hl2);
	  double s2 = (xyz22[1] - xyz21[1]) / (2.*hl2);
	  double c2 = (xyz22[2] - xyz21[2]) / (2.*hl2);
	  double dist2 = -xyz21[1] * c2 + xyz21[2] * s2;
	  double pitch2 = geom->WirePitch(0, 1, plane2, tpc);

	  assert(hitmap[tpc][plane1].size() <= hitmap[tpc][plane2].size());

	  // Loop over pairs of hits.

	  art::PtrVector<recob::Hit> hitvec;
	  hitvec.reserve(2);

	  for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator
		ihit1 = hitmap[tpc][plane1].begin();
	      ihit1 != hitmap[tpc][plane1].end(); ++ihit1) {

	    const art::Ptr<recob::Hit>& phit1 = ihit1->second;
	    unsigned short channel1 = phit1->Channel();

	    // Get endpoint coordinates of this wire.

	    unsigned int tpc1a, plane1a, wire1;
	    const geo::WireGeo& wgeo = geom->ChannelToWire(channel1, tpc1a, plane1a, wire1);
	    assert(tpc1a == tpc);
	    assert(plane1a == plane1);
	    double hl1 = wgeo.HalfL();
	    double xyz1[3];
	    double xyz2[3];
	    wgeo.GetCenter(xyz1, -hl1);
	    wgeo.GetCenter(xyz2, hl1);

	    // Find the plane2 wire numbers corresponding to the endpoints.

	    double wire21 = (-xyz1[1] * c2 + xyz1[2] * s2 - dist2) / pitch2;
	    double wire22 = (-xyz2[1] * c2 + xyz2[2] * s2 - dist2) / pitch2;

	    int wmin = std::min(wire21, wire22);
	    int wmax = std::max(wire21, wire22) + 1.;

	    for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator 
		  ihit2 = hitmap[tpc][plane2].lower_bound(wmin);
		ihit2 != hitmap[tpc][plane2].upper_bound(wmax); ++ihit2) {

	      const art::Ptr<recob::Hit>& phit2 = ihit2->second;

	      // Check current pair of hits for compatibility.
	      // By construction, hits should always have compatible views 
	      // and times, but may not have compatible mc information.

	      hitvec.clear();
	      hitvec.push_back(phit1);
	      hitvec.push_back(phit2);
	      bool ok = compatible(hitvec, timeOffset, maxDT, maxS);
	      if(ok) {

		// Add a space point.

		++n2;
		if(fFilter) {
		  sptkey_type key = &*phit2;
		  std::multimap<sptkey_type, SpacePointX>::iterator it = 
		    sptmap.insert(std::pair<sptkey_type, SpacePointX>(key, SpacePointX()));
		  sptkeys.insert(key);
		  SpacePointX& sptx = it->second;
		  sptx.goodness = fillSpacePoint(hitvec, timeOffset, sptx);
		}
		else {
		  spts.push_back(recob::SpacePoint());
		  fillSpacePoint(hitvec, timeOffset, spts.back());
		}
	      }
	    }
	  }
	}
      }
    }

    // If three-view space points are allowed, make a triple loop
    // over hits and produce space points for compatible triplets.

    if(nplane >= 3 && fMinViews <= 3) {

      // Loop over triplets of hits.

      art::PtrVector<recob::Hit> hitvec;
      hitvec.reserve(3);

      unsigned int plane1 = index[0];
      unsigned int plane2 = index[1];
      unsigned int plane3 = index[2];

      // Get angle, pitch, and offset of plane1 wires.

      const geo::WireGeo& wgeo1 = geom->Plane(plane1, tpc).Wire(0);
      double hl1 = wgeo1.HalfL();
      double xyz11[3];
      double xyz12[3];
      wgeo1.GetCenter(xyz11, -hl1);
      wgeo1.GetCenter(xyz12, hl1);
      double s1 = (xyz12[1] - xyz11[1]) / (2.*hl1);
      double c1 = (xyz12[2] - xyz11[2]) / (2.*hl1);
      double dist1 = -xyz11[1] * c1 + xyz11[2] * s1;
      double pitch1 = geom->WirePitch(0, 1, plane1, tpc);

      // Get angle, pitch, and offset of plane2 wires.

      const geo::WireGeo& wgeo2 = geom->Plane(plane2, tpc).Wire(0);
      double hl2 = wgeo2.HalfL();
      double xyz21[3];
      double xyz22[3];
      wgeo2.GetCenter(xyz21, -hl2);
      wgeo2.GetCenter(xyz22, hl2);
      double s2 = (xyz22[1] - xyz21[1]) / (2.*hl2);
      double c2 = (xyz22[2] - xyz21[2]) / (2.*hl2);
      double dist2 = -xyz21[1] * c2 + xyz21[2] * s2;
      double pitch2 = geom->WirePitch(0, 1, plane2, tpc);

      // Get angle, pitch, and offset of plane3 wires.

      const geo::WireGeo& wgeo3 = geom->Plane(plane3, tpc).Wire(0);
      double hl3 = wgeo3.HalfL();
      double xyz31[3];
      double xyz32[3];
      wgeo3.GetCenter(xyz31, -hl3);
      wgeo3.GetCenter(xyz32, hl3);
      double s3 = (xyz32[1] - xyz31[1]) / (2.*hl3);
      double c3 = (xyz32[2] - xyz31[2]) / (2.*hl3);
      double dist3 = -xyz31[1] * c3 + xyz31[2] * s3;
      double pitch3 = geom->WirePitch(0, 1, plane3, tpc);

      // Get sine of angle differences.

      double s12 = s1 * c2 - s2 * c1;   // sin(theta1 - theta2).
      double s23 = s2 * c3 - s3 * c2;   // sin(theta2 - theta3).
      double s31 = s3 * c1 - s1 * c3;   // sin(theta3 - theta1).

      // Loop over hits in plane1.

      for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator
	    ihit1 = hitmap[tpc][plane1].begin();
	  ihit1 != hitmap[tpc][plane1].end(); ++ihit1) {

	unsigned int wire1 = ihit1->first;
	const art::Ptr<recob::Hit>& phit1 = ihit1->second;
	unsigned short channel1 = phit1->Channel();

	// Get endpoint coordinates of this wire from plane1.

	unsigned int tpc1a, plane1a, wire1a;
	const geo::WireGeo& wgeo = geom->ChannelToWire(channel1, tpc1a, plane1a, wire1a);
	assert(tpc1a == tpc);
	assert(plane1a == plane1);
	assert(wire1a == wire1);
	double hl1 = wgeo.HalfL();
	double xyz1[3];
	double xyz2[3];
	wgeo.GetCenter(xyz1, -hl1);
	wgeo.GetCenter(xyz2, hl1);

	// Get corrected time and oblique coordinate of first hit.

	double t1 = phit1->PeakTime() - timeOffset[tpc][plane1];
	double u1 = wire1 * pitch1 + dist1;

	// Find the plane2 wire numbers corresponding to the endpoints.

	double wire21 = (-xyz1[1] * c2 + xyz1[2] * s2 - dist2) / pitch2;
	double wire22 = (-xyz2[1] * c2 + xyz2[2] * s2 - dist2) / pitch2;

	int wmin = std::min(wire21, wire22);
	int wmax = std::max(wire21, wire22) + 1.;

	for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator 
	      ihit2 = hitmap[tpc][plane2].lower_bound(wmin);
	    ihit2 != hitmap[tpc][plane2].upper_bound(wmax); ++ihit2) {

	  int wire2 = ihit2->first;
	  const art::Ptr<recob::Hit>& phit2 = ihit2->second;

	  // Get corrected time of second hit.

	  double t2 = phit2->PeakTime() - timeOffset[tpc][plane2];

	  // Check maximum time difference with first hit.

	  bool dt12ok = std::abs(t1-t2) <= maxDT;
	  if(dt12ok) {

	    // Test first two hits for compatibility before looping 
	    // over third hit.

	    hitvec.clear();
	    hitvec.push_back(phit1);
	    hitvec.push_back(phit2);
	    bool h12ok = compatible(hitvec, timeOffset, maxDT, maxS);
	    if(h12ok) {

	      // Get oblique coordinate of second hit.

	      double u2 = wire2 * pitch2 + dist2;

	      // Predict plane3 oblique coordinate and wire number.

	      double u3pred = (-u1*s23 - u2*s31) / s12;
	      double w3pred = (u3pred - dist3) / pitch3;
	      double w3delta = std::abs(maxS / (s12 * pitch3));
	      int w3min = std::ceil(w3pred - w3delta);
	      int w3max = std::floor(w3pred + w3delta);

	      for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator 
		    ihit3 = hitmap[tpc][plane3].lower_bound(w3min);
		  ihit3 != hitmap[tpc][plane3].upper_bound(w3max); ++ihit3) {

		int wire3 = ihit3->first;
		const art::Ptr<recob::Hit>& phit3 = ihit3->second;

		// Get corrected time of third hit.

		double t3 = phit3->PeakTime() - timeOffset[tpc][plane3];

		// Check time difference of third hit compared to first two hits.

		bool dt123ok = std::abs(t1-t3) <= maxDT && std::abs(t2-t3) <= maxDT;
		if(dt123ok) {

		  // Get oblique coordinate of third hit and check spatial separation.

		  double u3 = wire3 * pitch3 + dist3;
		  double S = s23 * u1 + s31 * u2 + s12 * u3;
		  bool sok = std::abs(S) <= maxS;
		  if(sok) {

		    // Test triplet for compatibility.

		    hitvec.clear();
		    hitvec.push_back(phit1);
		    hitvec.push_back(phit2);
		    hitvec.push_back(phit3);
		    bool h123ok = compatible(hitvec, timeOffset, maxDT, maxS);
		    if(h123ok) {

		      // Add a space point.

		      ++n3;
		      if(fFilter) {
			sptkey_type key = &*phit3;
			std::multimap<sptkey_type, SpacePointX>::iterator it = 
			  sptmap.insert(std::pair<sptkey_type, SpacePointX>(key, SpacePointX()));
			sptkeys.insert(key);
			SpacePointX& sptx = it->second;
			sptx.goodness = fillSpacePoint(hitvec, timeOffset, sptx);
		      }
		      else {
			spts.push_back(recob::SpacePoint());
			fillSpacePoint(hitvec, timeOffset, spts.back());
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // Do Filtering.

    if(fFilter) {

      // Transfer (some) space points from sptmap to spts.

      spts.reserve(spts.size() + sptkeys.size());

      // Loop over keys of space point map.
      // Space points that have the same key are candidates for filtering.

      for(std::set<sptkey_type>::const_iterator i = sptkeys.begin();
	  i != sptkeys.end(); ++i) {
	sptkey_type key = *i;

	// Loop over space points corresponding to the current key.
	// Choose the single best space point from among this group.

	double best_goodness = 0.;
	const SpacePointX* best_sptx = 0;

	for(std::multimap<sptkey_type, SpacePointX>::const_iterator j = sptmap.lower_bound(key);
	    j != sptmap.upper_bound(key); ++j) {
	  const SpacePointX& sptx = j->second;
	  if(best_sptx == 0 || sptx.goodness < best_goodness) {
	    best_sptx = &sptx;
	    best_goodness = sptx.goodness;
	  }
	}

	// Transfer best filtered space point to result vector.

	assert(best_sptx != 0);
	if(best_sptx != 0) {
	  spts.push_back(*best_sptx);
	  if(fMinViews <= 2)
	    ++n2filt;
	  else
	    ++n3filt;
	}
      }
    }
    else {
      n2filt = n2;
      n3filt = n3;
    }
  }

  debug << "\n2-hit space points = " << n2 << "\n"
	<< "3-hit space points = " << n3 << "\n"
	<< "2-hit filtered space points = " << n2filt << "\n"
	<< "3-hit filtered space points = " << n3filt;
}
