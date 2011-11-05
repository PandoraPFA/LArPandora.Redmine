////////////////////////////////////////////////////////////////////////////
// \version $Id: EndPoint2D.cxx,v 1.2 2010/06/19 22:20:12 spitz7 Exp $
//
// \brief Definition of vertex object for LArSoft
//
// \author spitz
//
////////////////////////////////////////////////////////////////////////////

#include "RecoBase/EndPoint2D.h"
#include <string>
#include <iostream>

namespace recob{

  //----------------------------------------------------------------------
  EndPoint2D::EndPoint2D() : 
    fDriftTime(-1.),
    fWireNum(-1),
    fID(-1),
    fStrength(-1.),
    fView(geo::kUnknown)
  {
  }

  //----------------------------------------------------------------------
  EndPoint2D::EndPoint2D(art::PtrVector<recob::Hit> &hits) :
    fDriftTime(-1.),
    fWireNum(-1),
    fID(-1),
    fStrength(-1.),
    fView(hits[0]->View()),
    fHits(hits)
  {
  }

  //----------------------------------------------------------------------
  EndPoint2D::~EndPoint2D()
  {
  }

  //----------------------------------------------------------------------
  ///if wire is < 0 then return all hits for a given plane and wire if specified
  ///just return those hits
 art::PtrVector<recob::Hit> EndPoint2D::Hits(unsigned int plane, unsigned int wire, bool allwires) const
  {
  
    unsigned int p = 0;
    unsigned int w = 0;
    unsigned int t = 0;

    art::PtrVector<recob::Hit> hits;
    art::ServiceHandle<geo::Geometry> geo;
    
    // \todo Do we need to check the TPC number here?
    for(size_t i = 0; i < fHits.size(); ++i){
      geo->ChannelToWire(fHits[i]->Wire()->RawDigit()->Channel(), t, p, w);
      if(p == plane && (w == wire || allwires) ) hits.push_back(fHits[i]);
    }
    
    return hits;
  }

  //----------------------------------------------------------------------
  art::PtrVector<recob::Hit> EndPoint2D::Hits(geo::View_t view) const
  {
    art::PtrVector<recob::Hit> hits;

    if(view == geo::kUnknown){
      std::cerr << "WARNING: view is geo::kUnknown, returning empty vector of hits" << std::endl;
      return hits;
    }

    for(size_t i = 0; i < fHits.size(); ++i){
      if(fHits[i]->View() == view) hits.push_back(fHits[i]);
    }// end loop over hits

    return hits;
  }

  //----------------------------------------------------------------------
  double EndPoint2D::Charge(geo::View_t view)
  {
    double charge = 0.;

    for(size_t i = 0; i < fHits.size(); ++i){
      art::Ptr<recob::Hit> hit = fHits[i];

      if(hit->View() == view){
	charge += hit->Charge();
      }///end if the hit is the correct view
    }///end loop over hits
  
    return charge;
  }

  //----------------------------------------------------------------------
  // ostream operator.  
  //

  std::ostream& operator<< (std::ostream& o, const EndPoint2D& ep)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << "Vertex2D ID "    << std::setw(5)  << std::right << ep.ID() 
      << " : View = "     << std::setw(3)  << std::right << ep.View()  
      << " Wire = "  << std::setw(7)  << std::right << ep.WireNum() 
      << " Time = "  << std::setw(9)  << std::right << ep.DriftTime();
      
    return o;
  }


}
