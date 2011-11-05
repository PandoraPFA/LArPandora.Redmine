////////////////////////////////////////////////////////////////////////
///
/// \file   SpacePointService.h
///
/// \brief  Service for generating space points from hits.
///
/// \author H. Greenlee 
///
/// This service calculates space points (recob::SpacePoint) from an 
/// unsorted collection of hits (recob::Hit).  The resulting space
/// points will contain one hit from from two or three views.
///
/// FCL parameters:
///
/// UseMC - Use MC truth information to find only true space points.
/// MaxDT - The maximum time difference (ticks) between any pair of hits.
/// MaxS  - The maximum 3-view wire separation parameter S (cm).
/// TimeOffsetU - Plane U time offset (ticks).
/// TimeOffsetV - Plane V time offset (ticks).
/// TimeOffsetW - Plane W time offset (ticks).
/// MinViews - Minimum number of views to make a space point (2 or 3).
/// EnableU - Use U view hits.
/// EnableV - Use V view hits.
/// EnableW - Use W view hits.
/// Filter - Filter space points flag.
///
/// The parameters fMaxDT and fMaxS are used to implement a notion of whether
/// the input hits are compatible with being a space point.  Parameter
/// MaxS is a cut on the 3-plane wire separation parameter S, which is
/// defined as follows:
///
/// S = sin(theta_vw)*u + sin(theta_wu)*v + sin(theta_uv)*w
///
/// where wire coordinates (u,v,w) are measured in cm with respect to a 
/// common origin.
///
/// The time offsets are subtracted from times embedded in hits before
/// comparing times in different planes, and before converting time
/// to distance. 
///
/// If enabled, filtering eliminates multiple space points with similar
/// times on the same wire of the most popluated plane.
///
/// There should eventually be a better way to specify time offsets.
///
////////////////////////////////////////////////////////////////////////

#ifndef SPACEPOINTSERVICE_H
#define SPACEPOINTSERVICE_H

#include <vector>
#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "RecoBase/SpacePoint.h"
#include "Simulation/SimChannel.h"

class TH1F;
namespace sim {
  class IDE;
}

namespace trkf {

  class SpacePointService {
  public:

    // Constructor.
    SpacePointService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

    // Destructor.
    ~SpacePointService();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    // Print constants obtained from geometry and properties services.
    void update() const;

    // Calculate time offsets.
    // Results stored in nested vector indexed by [tpc][plane]
    void fillTimeOffset(std::vector<std::vector<double> >& timeOffset) const;

    // Corrected time accessors.
    double correctedTime(const recob::Hit& hit,
			 const std::vector<std::vector<double> >& timeOffset) const;

    // Spatial separation of hits (zero if two or fewer).
    double separation(const art::PtrVector<recob::Hit>& hits) const;

    // Test whether the specified hits are compatible with a space point.
    // The last two arguments can be used to override the default cuts.
    bool compatible(const art::PtrVector<recob::Hit>& hits,
		    const std::vector<std::vector<double> >& timeOffset,
		    double maxDT = 0.,
		    double maxS = 0.) const;

    // Fill a single space point using the specified hits.
    // Hits are assumed to be compatible.
    // Returned value is a time goodness of fit (smaller is better, like chisquare).
    double fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
			  const std::vector<std::vector<double> >& timeOffset,
			  recob::SpacePoint& spt) const;

    // Fill a vector of space points from an unsorted collection of hits.
    // Space points are generated for all compatible combinations of hits.
    // The last two arguments can be used to override the default cuts.
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts,
			 double maxDT = 0.,
			 double maxS = 0.) const;

    // The following versions of the makeSpacePoints methods can be used
    // to fill only those space points that are compatible with mc truth
    // information contained in SimChannels.
    // The last two arguments can be used to override the default cuts.
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts,
			 const std::vector<const sim::SimChannel*>& simchans,
			 double maxDT = 0.,
			 double maxS = 0.) const;

  private:

    // Extended SpacePoint struct.
    // Contains extra information used for filtering (private type).

    struct SpacePointX : public recob::SpacePoint {
      SpacePointX() : goodness(0) {}
      double goodness;
    };

    // Configuration paremeters.

    bool fUseMC;            ///< Use MC truth information.
    double fMaxDT;          ///< Maximum time difference between planes.
    double fMaxS;           ///< Maximum space separation between wires.
    double fTimeOffsetU;    ///< Time offset corrections (U)
    double fTimeOffsetV;    ///< Time offset corrections (V)
    double fTimeOffsetW;    ///< Time offset corrections (W)
    int fMinViews;          ///< Mininum number of views per space point.
    bool fEnableU;          ///< Enable flag (U).
    bool fEnableV;          ///< Enable flag (V).
    bool fEnableW;          ///< Enable flag (W).
    bool fFilter;           ///< Filter flag.

    // Temporary variables.

    struct HitMCInfo
    {
      std::vector<int> trackIDs;       ///< Parent trackIDs.
      std::vector<double> xyz;         ///< Location of ionization (all tracks).
      std::vector<const recob::Hit*> pchit;   ///< Pointer to nearest neighbor hit (indexed by plane).
      std::vector<double> dist2;              ///< Distance to nearest neighbor hit (indexed by plane).
    };
    mutable std::map<const recob::Hit*, HitMCInfo> fHitMCMap;
  };
}

#endif
