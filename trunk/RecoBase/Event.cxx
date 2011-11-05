////////////////////////////////////////////////////////////////////////////
// \version $Id: Vertex.cxx,v 1.4 2010/06/10 16:21:31 antonm Exp $
//
// \brief Definition of event object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoBase/Event.h"
#include "TMath.h"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>

namespace recob{

  //----------------------------------------------------------------------
  Event::Event()
  {
  }

  //----------------------------------------------------------------------
  Event::Event(art::PtrVector<recob::Vertex> &vertices,
	       int id)
    : fVertices(vertices)
    , fID(id)
  {
    // sort the vertices, set the first one in z to be the primary
    fVertices.sort();
    fPrimary = fVertices[0];
  }

  //----------------------------------------------------------------------
  Event::~Event()
  {
  }

  //----------------------------------------------------------------------
  art::PtrVector<recob::Hit> Event::Hits() const
  {
    art::PtrVector<recob::Hit> hits;

    // loop over all the vertices and get their hits
    for(size_t v = 0; v < fVertices.size(); ++v){
      art::PtrVector<recob::Hit> hs = fVertices[v]->Hits();
      for(size_t h = 0; h < hs.size(); ++h) hits.push_back(hs[h]);
    }

    return hits;
  }

  //----------------------------------------------------------------------
  double Event::Energy() const
  {
    // loop over all vertex objects and get the 
    mf::LogWarning("Event") << "Event::Energy() is not yet defined.  Need to decide "
			    << " how to calculate energy of Prong/Vertex"
			    << " Return -999. for now.";
    
    return -999.;
  }

  //----------------------------------------------------------------------
  double Event::SigmaEnergy() const
  {
    // loop over all vertex objects and get the 
    mf::LogWarning("Event") << "Event::SigmaEnergy() is not yet defined.  Need to decide "
			    << " how to calculate uncertainty in energy of Prong/Vertex"
			    << " Return -999. for now.";
    
    return -999.;
  }

   //----------------------------------------------------------------------
  // Print function...to be called by << operator.  Facilitates overriding
  // by inheriting Track/Shower classes.
  //
   std::ostream& Event::Print(std::ostream& stream) const
   {
     stream << std::setprecision(5);
     stream << "Event " << this->fID << std::setw(5)  
	    << " #Vertices = " << std::setw(5) << std::right <<  this->fVertices.size()
	    << " Energy = " << this->Energy() << " +/-" << this->SigmaEnergy();
     
      return stream;
     
   }

  //----------------------------------------------------------------------
  // ostream operator.  
  //
  std::ostream& operator<< (std::ostream& o, const Event & a)
  {
     return a.Print(o);
  }


  //----------------------------------------------------------------------
  // < operator.  
  //
  bool operator < (const Event & a, const Event & b)
  {
    double xyza[3] = {0.};
    double xyzb[3] = {0.};
    a.PrimaryVertex()->XYZ(xyza);
    b.PrimaryVertex()->XYZ(xyzb);

    return xyza[2] < xyzb[2];

  }


}
