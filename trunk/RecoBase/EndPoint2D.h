////////////////////////////////////////////////////////////////////////////
// \version $Id: Vertex.h,v 1.2 2010/06/19 22:20:12 spitz7 Exp $
//
// \brief Definition of vertex object for LArSoft
//
// \author joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////////

#ifndef RECOB_ENDPOINT2D_H
#define RECOB_ENDPOINT2D_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "art/Persistency/Common/PtrVector.h"

#include "RecoBase/Hit.h"
#include "Geometry/geo.h"

namespace recob {
  
  class EndPoint2D  {

  public:
    
    EndPoint2D();  ///Default constructor
    explicit EndPoint2D(art::PtrVector<recob::Hit> &hits);
    ~EndPoint2D();
    void                           SetDriftTime(double drifttime) {fDriftTime = drifttime;}
    void                           SetWireNum(int wire)           {fWireNum = wire;       }
    void                           SetID(int ID)                  {fID=ID;                }
    void                           SetStrength(double Strength)   {fStrength=Strength;    }
    void                           SetView(geo::View_t view)      {fView = view;}
    art::PtrVector<recob::Hit> Hits(unsigned int plane, unsigned int wire=0, bool allwires=true) const;
    art::PtrVector<recob::Hit> Hits(geo::View_t view=geo::kUnknown)                              const;
    double                         Charge(geo::View_t view=geo::kUnknown);
    geo::View_t                    View()      const { return fView;      }
    double                         DriftTime() const { return fDriftTime; }     
    int                            WireNum()   const { return fWireNum;   }    
    int                            ID()        const { return fID;        }
    double                         Strength()  const { return fStrength;  }

    friend std::ostream& operator << (std::ostream& o, const EndPoint2D& c);
  
   private:
    double                     fDriftTime; // vertex's drift time
    int                        fWireNum;   // vertex's wire number
    int                        fID;        // vertex's ID
    double                     fStrength;  // vertex's strength   
    geo::View_t                fView;      // view for this vertex, geo::3D if 3D
    art::PtrVector<recob::Hit> fHits;      // vector of hits indices

  };
}

#endif //RECOB_ENDPOINT2D_H
