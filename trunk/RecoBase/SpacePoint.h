////////////////////////////////////////////////////////////////////////////
//
// Definition of SpacePoint class for LArSoft
//
// SpacePoints are 3D objects that contain pointers to Hits from multiple
// wireplanes that have been identified as matching.
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////////

#ifndef SPACEPOINT_H
#define SPACEPOINT_H

#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "art/Persistency/Common/PtrVector.h"

#include "RecoBase/Hit.h"
#include "Geometry/geo.h"

namespace recob {
  
  class SpacePoint {

  public:
    
    SpacePoint();  ///Default constructor 
    explicit SpacePoint(const art::PtrVector<recob::Hit> &hits); 
    ~SpacePoint();

    void SetID (int ID)           { fID=ID;}
    void SetXYZ(double *xyz)      { fXYZ[0] = xyz[0]; fXYZ[1] = xyz[1]; fXYZ[2] = xyz[2];}
      

    art::PtrVector<recob::Hit> Hits(unsigned int plane, 
                                    bool allhits=false)               const;
    art::PtrVector<recob::Hit> Hits(geo::View_t view=geo::kUnknown, 
                                    bool allhits=false)               const;
    int                        ID()                                   const { return fID; }  
    const double*              XYZ()                                  const { return fXYZ;}
    void                       PrintHits()                            const;  
     
    friend std::ostream& operator << (std::ostream& o, const SpacePoint & a);
    friend bool          operator <  (const SpacePoint & a, const SpacePoint & b);
    
  
  private:
    int                        fID;        ///SpacePoint ID
    double                     fXYZ[3];    ///position of SpacePoint in xyz
    art::PtrVector<recob::Hit> fHits;      ///ptrvector of hits in this SpacePoint

  };
}

#endif //SPACEPOINT_H
