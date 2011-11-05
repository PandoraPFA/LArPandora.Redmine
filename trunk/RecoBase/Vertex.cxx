////////////////////////////////////////////////////////////////////////////
// \version $Id: Vertex.cxx,v 1.4 2010/06/10 16:21:31 antonm Exp $
//
// \brief Definition of vertex object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoBase/Vertex.h"
#include "TMath.h"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>

namespace recob{

  bool sort_hits(art::Ptr<recob::Hit> const& a, art::Ptr<recob::Hit> const& b)
  {
    return a.get() < b.get();
  }

  bool duplicate_hits(art::Ptr<recob::Hit> const& a, art::Ptr<recob::Hit> const& b)
  {
    if(a->Channel()    != b->Channel()  ) return false;
    if(a->StartTime()  != b->StartTime()) return false;
    if(a->EndTime()    != b->EndTime()  ) return false;

    return true;
  }

  //----------------------------------------------------------------------
  Vertex::Vertex()
  {
  }

  //----------------------------------------------------------------------
  Vertex::Vertex(art::PtrVector<recob::Track>  &tracks, 
		 art::PtrVector<recob::Shower> &showers, 
		 double *xyz,
		 int id) 
    : fTracks (tracks)
    , fShowers(showers)
    , fID(id)
  {
    fXYZ[0] = xyz[0];
    fXYZ[1] = xyz[1];
    fXYZ[2] = xyz[2];
  }

  //----------------------------------------------------------------------
  Vertex::~Vertex()
  {
  }

  //----------------------------------------------------------------------
  void Vertex::XYZ(double *xyz) const
  {
    xyz[0] = fXYZ[0];
    xyz[1] = fXYZ[1];
    xyz[2] = fXYZ[2];
    
    return;
  }

  //----------------------------------------------------------------------
  art::PtrVector<recob::Hit> Vertex::Hits() const
  {
    // make an std::vector of the hits that we can sort
    // and drop duplicate entries from
    std::vector< art::Ptr<recob::Hit> > hits;

    // loop over all the showers and tracks and get their hits
    for(size_t t = 0; t < fTracks.size(); ++t){
      art::PtrVector<recob::Hit> hs = fTracks[t]->Hits();
      for(size_t h = 0; h < hs.size(); ++h) hits.push_back(hs[h]);
    }
    for(size_t s = 0; s < fShowers.size(); ++s){
      art::PtrVector<recob::Hit> hs = fShowers[s]->Hits();
      for(size_t h = 0; h < hs.size(); ++h) hits.push_back(hs[h]);
    }

    std::sort(hits.begin(), hits.end(), sort_hits);
    std::vector< art::Ptr<recob::Hit> >::iterator it = std::unique(hits.begin(), hits.end(), duplicate_hits);
    hits.resize(it - hits.begin());

    // now fill an art::PtrVector to return
    art::PtrVector<recob::Hit> hitvec;
    for(size_t h = 0; h < hits.size(); ++h) hitvec.push_back(hits[h]);

    return hitvec;
  }

  //----------------------------------------------------------------------
  // Print function...to be called by << operator.  Facilitates overriding
  // by inheriting Track/Shower classes.
  //
   std::ostream& Vertex::Print(std::ostream& stream) const
   {
     stream << std::setprecision(5);
     stream << "Vertex ID    "  << this->fID << std::setw(5)  
	    << " : #Tracks  = " << std::setw(5) << std::right <<  this->fTracks.size()
	    << " #Showers = "   << std::setw(5) << std::right <<  this->fShowers.size()
	    << " (x,y,z)  = ("  << this->fXYZ[0] << "," << this->fXYZ[1] << "," 
	    << this->fXYZ[2] << ")";
     
      return stream;
     
   }

  //----------------------------------------------------------------------
  // ostream operator.  
  //
  std::ostream& operator<< (std::ostream& o, const Vertex & a)
  {
     return a.Print(o);
  }


  //----------------------------------------------------------------------
  // < operator.  
  //
  bool operator < (const Vertex & a, const Vertex & b)
  {
    double xyza[3] = {0.};
    double xyzb[3] = {0.};
    a.XYZ(xyza);
    b.XYZ(xyzb);

    return xyza[2] < xyzb[2];

  }


}
