////////////////////////////////////////////////////////////////////////////
// \version $Id: Prong.h,v 1.4 2010/06/10 16:21:31 antonm Exp $
//
// \brief Definition of prong object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#ifndef PRONG_H
#define PRONG_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "Geometry/geo.h"
#include "art/Persistency/Common/PtrVector.h"


namespace recob {
  
    class Prong  {

    public:
        
        Prong();  // Default constructor
        explicit Prong(art::PtrVector<recob::Cluster> &clusters, 
                       std::vector<recob::SpacePoint> spacepoints = std::vector<recob::SpacePoint>());// default initialize SpacePoint vector
        ~Prong();
        
        void                              SetID(int ID){fID=ID;}
        void                              SetDirection(double *dcosStart, double *dcosEnd);

	art::PtrVector<recob::Hit>            Hits(geo::View_t=geo::kUnknown)            const;
        art::PtrVector<recob::Cluster>        Clusters(geo::View_t view=geo::kUnknown)   const;
        const std::vector<recob::SpacePoint>& SpacePoints()                              const  {return fSpacePoints;}
        
        void                              Extent(std::vector<double> &xyzStart, std::vector<double> &xyzEnd) const; 
        void                              Direction(double *dcosStart, double *dcosEnd) const;  
        double                            ProjectedLength(geo::View_t view)             const;
        int                               ID()                                          const  {return fID;}  
        double                            MaxTransverseWidth(geo::Coord_t coordinate)   const; 
        double                            dSMaxTransverseWidth()                        const  {return fDistanceMaxWidth;}
        double                            Theta() const;
        double                            Phi()   const;
  
     
   
        friend bool          operator <   (const Prong & a, const Prong & b);
        
    protected:
        virtual std::ostream&          Print(std::ostream& stream) const;
        friend  std::ostream& operator <<   (std::ostream& o, const Prong & a);
        
        art::PtrVector<recob::Cluster> fClusters;       ///< collection of clusters
        std::vector<recob::SpacePoint> fSpacePoints;    ///< vector of spacepoints
        
    private:
        
        int    fID;                    ///< prong's ID
        double fDCosStart[3];          ///< direction cosines at start of prong
        double fSigmaDCosStart[3];     ///< uncertainting on initial direction cosines
        double fDCosEnd[3];            ///< direction cosines at end of prong  
        double fMaxTransverseWidth[2]; ///< maximum width of the prong in the x(0) and y(0) directions
        double fDistanceMaxWidth;      ///< distance from the start of the prong to its maximum width
        
    };
}

#endif // PRONG_H
