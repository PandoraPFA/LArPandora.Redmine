////////////////////////////////////////////////////////////////////////
/// \file  TPCGeo.h
/// \brief Encapsulate the construction of a single detector plane
///
/// \version $Id: TPCGeo.h,v 1.7 2009/12/01 21:07:51 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_TPCGEO_H
#define GEO_TPCGEO_H
#include <vector>
#include "TGeoVolume.h"
class TGeoNode;
class TGeoHMatrix;

namespace geo {
  class PlaneGeo;

  typedef enum driftdir {
    kUnknownDrift, ///< drift direction is unknown
    kPosX,         ///< drift towards positive X values			
    kNegX, 	   ///< drift towards negative X values			
    kPosY 	   ///< drift towards positive Y values (Bo only for now)
  } DriftDirection_t;

  //......................................................................
  /// Geometry information for a single tpc
  class TPCGeo {
  public:
    // Construct a representation of a single plane of the detector
    TPCGeo(std::vector<const TGeoNode*>& path, int depth);
    ~TPCGeo();

    // Number of planes in this tpc
    unsigned int Nplanes()                                      const { return fPlanes.size();   }       
								                                         
    // Return the iplane'th plane in the tpc. 			                                         
    const PlaneGeo& Plane(unsigned int iplane)                	const;				       
    								                                         
    double            ActiveHalfWidth()                         const; // half width of the active TPC 	       
    double            ActiveHalfHeight()          		const; // half height of the active TPC	       
    double            ActiveLength()              		const; // length of the active TPC     	  
    double            HalfWidth()                             	const; // half width of the total TPC 	       
    double            HalfHeight()          		      	const; // half height of the total TPC	       
    double            Length()              		      	const; // length of the total TPC     	  
    double            ActiveMass()                            	const { return fActiveVolume->Weight(); }
    const TGeoVolume* ActiveVolume()                          	const { return fActiveVolume; }          
    const TGeoVolume* TotalVolume()                          	const { return fTotalVolume; }          

    DriftDirection_t  DriftDirection()                          const { return fDriftDirection; }


    unsigned int      Nchannels()                               const { return fChannelMap.size(); }
    unsigned int      NearestChannel(double* worldPos)          const;
    unsigned int      PlaneWireToTPCChannel(unsigned int plane,
					    unsigned int wire)  const; // convert plane, wire to channel in TPC
    const WireGeo&    ChannelToWire(unsigned int  channel,
				    unsigned int &plane, 
				    unsigned int &wire)         const; // convert tpc channel # to plane, wire
    const double*     PlaneLocation(unsigned int p)             const; 
    double            Plane0Pitch(unsigned int p)               const;
    double            PlanePitch(unsigned int p1=0,
				 unsigned int p2=1)             const;
    double            WirePitch(unsigned int w1=0,
				unsigned int w2=1,
				unsigned int p=0)               const;

    // Transform point from local plane frame to world frame
    void LocalToWorld(const double* tpc, double* world)         const;
    
    // Transform direction vector from local to world
    void LocalToWorldVect(const double* tpc, double* world)     const;

    // Again, with TVectors
    const TVector3 LocalToWorld( const TVector3& local )        const;

    // Transform point from world frame to local tpc frame
    void WorldToLocal(const double* world, double* tpc)         const;

    // Transform direction vector from world to local
    void WorldToLocalVect(const double* world, double* tpc)     const;
    
    // Again, with TVectors
    const TVector3 WorldToLocal( const TVector3& world )        const;
    
  private:
    
    void FindPlane(std::vector<const TGeoNode*>& path,
		   unsigned int depth);
    void MakePlane(std::vector<const TGeoNode*>& path, 
		   int depth);
    
  private:

    typedef std::pair<unsigned int, unsigned int>   PlaneWirePair;
    typedef std::map<unsigned short, PlaneWirePair> TPCChannelMap;

    TGeoHMatrix*                       fGeoMatrix;      ///< TPC to world transform				      	
    std::vector<PlaneGeo*> 	       fPlanes;         ///< List of planes in this plane			      	
    TGeoVolume*            	       fActiveVolume;   ///< Active volume of LAr, called volTPCActive in GDML file 
    TGeoVolume*                        fTotalVolume;    ///< Total volume of TPC, called volTPC in GDML file	
    TPCChannelMap          	       fChannelMap;     ///< Map channel number in TPC to plane and wire numbers    
    DriftDirection_t       	       fDriftDirection; ///< Direction of the electron drift in the TPC		
    std::vector<double>    	       fPlane0Pitch;    ///< Pitch between planes                                   
    std::vector< std::vector<double> > fPlaneLocation;  ///< xyz locations of planes in the TPC

  };
}

#endif
////////////////////////////////////////////////////////////////////////
