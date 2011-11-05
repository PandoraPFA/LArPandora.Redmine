////////////////////////////////////////////////////////////////////////////
// \version $Id: Cluster.h,v 1.6 2010/06/12 21:46:34 spitz7 Exp $
//
// \brief Definition of cluster object for LArSoft
//
// \author mitchell.soderberg@yale.edu
//
////////////////////////////////////////////////////////////////////////////

#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "math.h"

#include "art/Persistency/Common/PtrVector.h"

#include "RecoBase/Hit.h"
#include "Geometry/geo.h"

namespace recob {
  
  class Cluster {

  public:
    
    Cluster();  ///Default constructor
    explicit Cluster(art::PtrVector<recob::Hit> &hits,
		     double startWire, double sigmaStartWire,
		     double startTime, double sigmaStartTime,
		     double endWire, double sigmaEndWire, 
		     double endTime, double sigmaEndTime,
		     double dTdW, double sigmadTdW,
		     double dQdW, double sigmadQdW,
		     int id); 
    ~Cluster();


    art::PtrVector<recob::Hit> Hits(unsigned int wire=0, bool allwires=true) const;
    double                Charge()        const;
    geo::View_t           View()          const { return fHits[0]->View(); }
    double                dTdW()          const { return fdTdW;    	   }
    double                dQdW()          const { return fdQdW;    	   }
    double                SigmadTdW()     const { return fSigmadTdW;       }
    double                SigmadQdW()     const { return fSigmadQdW;       }
    std::vector<double>   StartPos()      const { return fStartPos;	   }
    std::vector<double>   EndPos()    	  const { return fEndPos;  	   }
    std::vector<double>   SigmaStartPos() const { return fSigmaStartPos;   }
    std::vector<double>   SigmaEndPos()   const { return fSigmaEndPos;     }
    int                   ID()            const { return fID;              }
    void                  PrintHits();  

    Cluster              operator +  (Cluster);
    friend std::ostream& operator << (std::ostream& o, const Cluster& c);
    friend bool          operator <  (const Cluster & a, const Cluster & b);
    
  
   private:
    double                     fdTdW;           ///< slope of cluster in tdc vs wire	     
    double                     fdQdW;      	///< slope of cluster in charge vs wire   
    double                     fSigmadTdW;      ///< slope of cluster in tdc vs wire	     
    double                     fSigmadQdW;      ///< slope of cluster in charge vs wire   
    std::vector<double>        fStartPos;  	///< start of cluster in (wire, tdc) plane
    std::vector<double>        fEndPos;    	///< start of cluster in (wire, tdc) plane
    std::vector<double>        fSigmaStartPos;  ///< start of cluster in (wire, tdc) plane
    std::vector<double>        fSigmaEndPos;    ///< start of cluster in (wire, tdc) plane
    int                        fID;             ///< cluster's ID
    art::PtrVector<recob::Hit> fHits;           ///< vector of hits in this cluster

  };
}

#endif //CLUSTER_H
