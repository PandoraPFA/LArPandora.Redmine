////////////////////////////////////////////////////////////////////////////
//
// Implementation of SpacePoint class for LArSoft
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////////

#include "RecoBase/SpacePoint.h"

namespace recob{

  //----------------------------------------------------------------------
  SpacePoint::SpacePoint() : 
    fID(-1)
  {
  }

  //----------------------------------------------------------------------
  SpacePoint::SpacePoint(const art::PtrVector<recob::Hit> &hits) :
    fID(-1),
    fHits(hits)
  {
  }

  //----------------------------------------------------------------------
  SpacePoint::~SpacePoint()
  {
  }

  //----------------------------------------------------------------------
  ///just return PtrVector of Hits associated with a SpacePoint from requested plane.
  // if "allhits" is provided as true, return all Hits associated with SpacePoint.
  // \todo make this method also require a tpc to be passed to it
  art::PtrVector<recob::Hit> SpacePoint::Hits(unsigned int plane, bool allhits) const
  {
    art::PtrVector<recob::Hit> hits;
    art::ServiceHandle<geo::Geometry> geo;

    unsigned int p = 0;
    unsigned int w = 0;
    unsigned int t = 0;

    for(size_t i = 0; i < fHits.size(); ++i){
      geo->ChannelToWire(fHits[i]->Channel(), t, p, w);
      if(p == plane || allhits)  hits.push_back(fHits[i]);
    }
    
    return hits;
  }

  //----------------------------------------------------------------------
  ///just return PtrVector of Hits associated with a SpacePoint from requested View.
  art::PtrVector<recob::Hit> SpacePoint::Hits(geo::View_t view, bool allhits) const
  {
    art::PtrVector<recob::Hit> hits;

    if(view == geo::kUnknown){
      std::cerr << "WARNING: view is geo::kUnknown, returning empty vector of hits" 
		<< std::endl;
      return hits;
    }

    for(size_t i = 0; i < fHits.size(); ++i){
      if(fHits[i]->View() == view || allhits) hits.push_back(fHits[i]);
    }

    return hits;
  }

  //----------------------------------------------------------------------
  // Print information for all Hits in a SpacePoints.  
  //
  void SpacePoint::PrintHits() const
  {
    for(size_t i = 0; i < fHits.size(); ++i) 
      std::cout << "Hit #" << i << " : " << *fHits[i] << std::endl;
    return;
  }
  
  //----------------------------------------------------------------------
  // ostream operator.  
  //

  std::ostream& operator<< (std::ostream& o, const SpacePoint & a)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << " SpacePoint ID " << std::setw(5) << std::right << a.ID() 
      << " #Hits = "       << std::setw(5) << std::right << a.Hits(-1).size()
      << " (X,Y,Z) = ("    << std::setw(5) << std::right << a.XYZ()[0]
      << " , "             << std::setw(5) << std::right << a.XYZ()[1]
      << " , "             << std::setw(5) << std::right << a.XYZ()[2]
      << ")" ;

    return o;
  }


  //----------------------------------------------------------------------
  // < operator.  
  //
  bool operator < (const SpacePoint & a, const SpacePoint & b)
  {
    if(a.ID() != b. ID())
      return a.ID()<b.ID();

    return false; //They are equal

  }

}

