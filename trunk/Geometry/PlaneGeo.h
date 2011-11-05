////////////////////////////////////////////////////////////////////////
/// \file  PlaneGeo.h
/// \brief Encapsulate the construction of a single detector plane
///
/// \version $Id: PlaneGeo.h,v 1.7 2009/12/01 21:07:51 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_PLANEGEO_H
#define GEO_PLANEGEO_H
#include <vector>
class TGeoNode;
class TGeoHMatrix;

namespace geo {
  class WireGeo;
  /// Enumerate the possible plane projections
  typedef enum _plane_proj {
    kUnknown, ///<unknown view
    kU,       ///< planes which measure U
    kV,       ///< planes which measure V
    kW,       ///< planes which measure W (third view for Bo, MicroBooNE, etc)
    k3D       ///< 3 dimensional objects, potentially hits, clusters, prongs, etc
  } View_t;
  typedef enum _plane_orient {
    kHorizontal, ///< planes that are in the horizontal plane (depricated as of 8/3/11 bjr)
    kVertical, ///< planes that are in the vertical plane (ie ArgoNeuT)
  } Orient_t;
  typedef enum _plane_sigtype {
    kMysteryType, ///< who knows?
    kInduction,   ///< signal from induction planes
    kCollection   ///< signal from collection planes
  } SigType_t;


  //......................................................................
  
  /// Geometry information for a single readout plane
  class PlaneGeo {
  public:
    /// Construct a representation of a single plane of the detector
    PlaneGeo(std::vector<const TGeoNode*>& path, int depth);
    ~PlaneGeo();

    /// Number of wires in this plane
    unsigned int Nwires()                                     const { return fWire.size();   }

    /// Return the iwire'th wire in the plane. 
    const WireGeo& Wire(unsigned int iwire)                   const { return *fWire[iwire];  }
    
    /// Which coordinate does this plane measure
    View_t View()                                             const { return fView;          }
    
    /// What is the orienation of the plane
    Orient_t Orientation()                                    const { return fOrientation;   }

    ///What is the signal type for the plane
    SigType_t SignalType()                                    const { return fSignalType;    }

    void SetSignalType(geo::SigType_t sigtype)                      { fSignalType = sigtype; }

    /// Transform point from local plane frame to world frame
    void LocalToWorld(const double* plane, double* world)     const;
    
    /// Transform direction vector from local to world
    void LocalToWorldVect(const double* plane, double* world) const;

    // Again, with TVectors
    const TVector3 LocalToWorld( const TVector3& local )      const;

    /// Transform point from world frame to local plane frame
    void WorldToLocal(const double* world, double* plane)     const;

    /// Transform direction vector from world to local
    void WorldToLocalVect(const double* world, double* plane) const;
    
    // Again, with TVectors
    const TVector3 WorldToLocal( const TVector3& world )      const;
    
  private:
    
    void FindWire(std::vector<const TGeoNode*>& path,
		  unsigned int depth);
    void MakeWire(std::vector<const TGeoNode*>& path, 
		  int depth);
    
  private:
    TGeoHMatrix*          fGeoMatrix;   ///< Plane to world transform
    View_t                fView;        ///< Does this plane measure U,V, or W?
    Orient_t              fOrientation; ///< Is the plane vertical or horizontal?
    SigType_t             fSignalType;  ///< Is the plane induction or collection?
    std::vector<WireGeo*> fWire;        ///< List of wires in this plane
  };
}

#endif
////////////////////////////////////////////////////////////////////////
