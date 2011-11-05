////////////////////////////////////////////////////////////////////////////
// \version $Id: Track.h,v 1.5 2010/02/15 20:32:46 brebel Exp $
//
// \brief Definition of track object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#ifndef TRACK_H
#define TRACK_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "RecoBase/Prong.h"

namespace recob {
  
  class Track :public Prong {

  public:
    
    Track();  ///Default constructor
    explicit Track(art::PtrVector<recob::Cluster> &clusters,
		   std::vector<recob::SpacePoint> spacepoints = std::vector<recob::SpacePoint>());//default initialize SpacePoint vector
    ~Track();
    
     std::ostream& Print(std::ostream& stream) const;

  };
}

#endif // TRACK_H
