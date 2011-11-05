////////////////////////////////////////////////////////////////////////////
// \version $Id: Prong.h,v 1.4 2010/06/10 16:21:31 antonm Exp $
//
// \brief Definition of vertex object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#ifndef RB_VERTEX_H
#define RB_VERTEX_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "RecoBase/Shower.h"
#include "RecoBase/Track.h"
#include "art/Persistency/Common/PtrVector.h"


namespace recob {
  
  class Vertex  {

  public:
    
    Vertex();  // Default constructor
    explicit Vertex(art::PtrVector<recob::Track>  &tracks,
		    art::PtrVector<recob::Shower> &showers,
		    double *xyz,
		    int id=-999);
    ~Vertex();

    art::PtrVector<recob::Hit>             Hits()           const;
    const art::PtrVector<recob::Track>&    Tracks()         const { return fTracks;  }
    const art::PtrVector<recob::Shower>&   Showers()        const { return fShowers; }
    void                                   XYZ(double *xyz) const;
    const int                              ID()             const { return fID;      }

    friend bool          operator <   (const Vertex & a, const Vertex & b);

  protected:
     virtual std::ostream&          Print(std::ostream& stream) const;
     friend  std::ostream& operator <<   (std::ostream& o, const Vertex & a);

    
  private:
    
    art::PtrVector<recob::Track>  fTracks;   ///< collection of Tracks
    art::PtrVector<recob::Shower> fShowers;  ///< collection of Showers
    double                        fXYZ[3];   ///< location of vertex
    int                           fID;       ///< id number for vertex
  };
}

#endif // RB_VERTEX_H
