////////////////////////////////////////////////////////////////////////
/// \file  Geometry.h
/// \brief Encapsulate the geometry of one entire detector
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
///
/// Revised <seligman@nevis.columbia.edu> 29-Jan-2009
///         Revise the class to make it into more of a general detector
///         interface.
///

#ifndef GEO_GEOMETRY_H
#define GEO_GEOMETRY_H

#include <TString.h>
#include <TVector3.h>
#include <Rtypes.h>

#include <vector>
#include <map>
#include <cstring>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

class TGeoManager;
class TGeoVolume;
class TGeoNode;
class TGeoMaterial;

namespace geo {

  // Foward declarations within namespace.
  class TPCGeo;
  class PlaneGeo;
  class WireGeo;

  typedef enum coordinates {
    kX, ///< X coordinate
    kY, ///< Y coordinate
    kZ  ///< Z coordinate
  } Coord_t;

  typedef enum detid {
    kBo,         ///< Bo id
    kArgoNeuT,   ///< ArgoNeuT id
    kMicroBooNE, ///< MicroBoone id
    kLBNE        ///< LBNE id
  } DetId_t;

  // The geometry of one entire detector.
  class Geometry 
  {
  public:

    // Access to the single instance of thie class.  For example:
    Geometry(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~Geometry();

    void preBeginRun(art::Run const& run);

    unsigned int Nchannels()                                      const;
    unsigned int NTPC()                                           const { return fTPCs.size();   }
    // Number of wire planes in the LAr TPC.
    unsigned int Nplanes(unsigned int tpc=0)                      const;
    // Number of wires in plane "i".
    unsigned int Nwires(unsigned int p, unsigned int tpc=0)       const;
    
    const TPCGeo&       TPC(unsigned int tpc)                     const;
    const TPCGeo&       PositionToTPC(double *worldLoc)           const; // return the TPCGeo object containing the world position worldLoc
    const PlaneGeo&     Plane(unsigned int p, unsigned int tpc=0) const;
    const WireGeo&      ChannelToWire(unsigned int  channel,
				      unsigned int &tpc,
				      unsigned int &plane, 
				      unsigned int &wire)         const; // convert channel # to tpc, plane, wire
    void                ChannelToWireFast(unsigned int  channel,
					  unsigned int &tpc,
					  unsigned int &plane, 
					  unsigned int &wire)     const; // convert channel # to tpc, plane, wire
                                                                         //  assuming heirachical numbering scheme
    unsigned int        PlaneWireToChannel(unsigned int plane, 
					   unsigned int wire,
					   unsigned int tpc=0)    const; // convert plane, wire to channel

    unsigned int        PlaneWireToChannelFast(unsigned int plane, 
					   unsigned int wire,
					   unsigned int tpc=0)    const; // convert plane, wire to channel
                                                                         //  assuming heirachical numbering scheme
    unsigned int        NearestChannel(double *worldLoc)          const; // find the nearest channel to 
                                                                         // input world coordinates
    unsigned int        NearestChannelFast(double * worldLoc, 
					   unsigned int PlaneNo, 
					   unsigned int TPCNo=0)  const; 
             

    const TGeoMaterial* Material(double x, 
				 double y, 
				 double z)                        const;				     
    double              DetHalfWidth(unsigned int tpc=0)          const; // half width of the TPC	     
    double              DetHalfHeight(unsigned int tpc=0)         const; // half height of the TPC     
    double              DetLength(unsigned int tpc=0)             const; // length of the TPC	     
    double              CryostatHalfWidth()                       const; // half width of the cryostat 
    double              CryostatHalfHeight()                      const; // half height of the cryostat
    double              CryostatLength()                          const; // length of the cryostat     
    void                CryostatBoundaries(double* boundaries)    const; // boundaries of cryostat, 3 pairs of +/- coord
    double              PlanePitch(unsigned int p1=0,                    // distance between planes    
				   unsigned int p2=1,		                                       
				   unsigned int tpc=0)            const; // p1 < p2		     
    double              WirePitch(unsigned int w1=0,                     // distance between wires     
				  unsigned int w2=1,                     // on the same plane	     
				  unsigned int plane=0,		                                       
				  unsigned int tpc=0)             const; // w1 < w2                    
    

    void                WorldBox(double* xlo, 
				 double* xhi,
				 double* ylo, 
				 double* yhi,
				 double* zlo, 
				 double* zhi)                     const; // volume box
    double              TotalMass(const char* vol="volCryostat")  const; // total mass of the 
                                                                         // specified volume
    double              MassBetweenPoints(double *p1, 
					  double *p2)             const; // mass between two points 
                                                                         // in the world

    // A typical y-position value at the surface (where earth meets air)
    // for this detector site
    //
    // \returns typical y position at surface in units of cm 
    double              SurfaceY()                                const { return fSurfaceY; }

    // Access to the ROOT geometry description.
    TGeoManager*        ROOTGeoManager()                          const;

    // The full directory path to the GDML file that was the source
    // of the detector geometry.
    std::string         GetGDMLPath()                             const { return fGDMLfile; }    
    std::string         GetROOTPath()                          	  const { return fROOTfile; }    
    std::string         ROOTFile()                                const { return fROOTfile; }
    std::string         GDMLFile()                                const { return fGDMLfile; }
    // The name of the detector.				                                 
    const TString       GetDetectorName()                      	  const { return fDetectorName; }

    // There are some issues that require detector-specific queries.
    // This method returns an enumerated type that can be tested in those cases
    geo::DetId_t DetId()                                       	  const {return fDetId; }


    // The Geant4 simulation needs to know the name of the world volume.
    const std::string GetWorldVolumeName()                        const;

    // The Geant4 simulation needs to know the name of the LAr TPC volume.
    const std::string GetLArTPCVolumeName(unsigned int tpc=0)     const;

    // The event display needs to know the name of the cryostat.
    const std::string GetCryostatVolumeName()                     const;

    // As of Aug-2009, the origin of the co-ordinate system for
    // ArgoNEUT and MicroBooNE is for z=0 and y=0 to be at the center
    // of the front face of the detector, but x=0 to be the edge of
    // the TPC.  This is convenient for read-out, but a pain for
    // simulation.  This method returns the center of the front face
    // of the TPC in the world co-ordinate system, making it easier
    // to write detector-independent simulation code.
    const TVector3 GetTPCFrontFaceCenter()                        const;

    // Name of the deepest volume containing the point xyz
    // returns volume containing the origin by default
    const std::string VolumeName(TVector3 point);

    // Name of the deepest material containing the point xyz
    // returns material of the origin by default
    const std::string MaterialName(TVector3 point);

    // The following functions are utilized to determine if two wires 
    // in the TPC intersect or not, and if they do then 
    // determine the coordinates of the intersection.
    // Starting point of wire is end with lower z-coordinate.
    bool ValueInRange(double value, 
		      double min, 
		      double max);
    void WireEndPoints(unsigned int tpc,
		       unsigned int plane, 
		       unsigned int wire, 
		       double *xyzStart, 
		       double *xyzEnd);
    bool ChannelsIntersect(unsigned short c1, 
			   unsigned short c2, 
			   double &y, 
			   double &z);
    void IntersectionPoint(unsigned int wire1, 
			   unsigned int wire2, 
			   unsigned int plane1, 
			   unsigned int plane2,
			   unsigned int tpc,
                           double start_w1[3], 
			   double end_w1[3], 
			   double start_w2[3], 
			   double end_w2[3], 
                           double &y, 
			   double &z);
    void IntersectionPoint(unsigned int wire1, 
			   unsigned int wire2, 
			   unsigned int plane1, 
			   unsigned int plane2,
			   unsigned int tpc,
                           double &y, 
			   double &z);
    bool CrossesVol(double xyz[], 
		    double dxyz[], 
		    double point[]); 


  private:
    
    void LoadGeometryFile(std::string gdmlfile, 
			  std::string rootfile);

    void FindTPC(std::vector<const TGeoNode*>& path,
		 unsigned int depth);
    void MakeTPC(std::vector<const TGeoNode*>& path, int depth);

    //types to handle the mapping of plane/wire to channel and back again
    typedef std::pair<unsigned int,  unsigned int>  TPCPlaneWireMap; 
    typedef std::map<unsigned int, TPCPlaneWireMap> ChannelMap;

    std::vector<TPCGeo*>   fTPCs;             ///< The detector planes
    double                 fSurfaceY;         ///< The point where air meets earth for this detector.
    TString                fDetectorName;     ///< Name of the detector.
    std::string            fGDMLfile;         ///< The GDML file used for the detector geometry - full path for now
    std::string            fROOTfile;         ///< The GDML file used for the detector geometry - full path for now
    double                 fMinWireZDist;     ///< Minimum distance in Z from a point in which to look for the closest wire
    bool                   fDisableWiresInG4; ///< If set true, supply G4 with GDMLfileNoWires rather than GDMLfile
    bool                   fForceUseFCLOnly;  ///< Force Geometry to only use the geometry files specified in the fcl file 
    geo::DetId_t           fDetId;            ///< Detector type
    ChannelMap             fChannelMap;       ///< map from channel # to a plane wire pair
    int                    fNBadVol;
  };

} // namespace geo

#endif // GEO_GEOMETRY_H
