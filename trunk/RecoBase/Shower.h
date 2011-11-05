////////////////////////////////////////////////////////////////////////////
// \version $Id: Shower.h,v 1.2 2010/02/15 20:32:46 brebel Exp $
//
// \brief Definition of shower object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#ifndef SHOWER_H
#define SHOWER_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "RecoBase/Prong.h"

namespace recob {
  
  class Shower : public Prong{

  public:
    
    Shower();  ///Default constructor
    explicit Shower(art::PtrVector<recob::Cluster> &clusters,
		    std::vector<recob::SpacePoint> spacepoints = std::vector<recob::SpacePoint>());//default initialize SpacePoint vector
    ~Shower();

    private:
            
  };
}

#endif // SHOWER_H
