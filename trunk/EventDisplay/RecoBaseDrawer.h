/// \file    RecoBaseDrawer.h
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  messier@indiana.edu
/// \version $Id: RecoBaseDrawer.h,v 1.3 2010/11/11 22:47:20 p-novaart Exp $
#ifndef EVD_RECOBASEDRAWER_H
#define EVD_RECOVASEDRAWER_H
#include <vector>
#ifdef __CINT__
namespace art { 
  class Event;
  class PtrVector;
  class Ptr;
}
#else
#include "art/Persistency/Common/PtrVector.h"
#endif

class TH1F;
namespace evdb     { class View2D;   }
namespace evdb     { class View3D;   }
namespace geo      { class Geometry; }
namespace recob { 
  class Hit;     
  class Cluster;
  class Prong;       
  class Track;       
  class Shower;      
  class Wire;
  class Vertex;
  class Event;
}

namespace evd {
  /// Aid in the rendering of RecoBase objects
  class RecoBaseDrawer {
  public:
    RecoBaseDrawer();
    ~RecoBaseDrawer();

  public:

    void Wire2D(const art::Event& evt,
		evdb::View2D*     view,
		unsigned int      plane);
    void Hit2D(const art::Event& evt,
	       evdb::View2D*     view,
	       unsigned int      plane);
    void Hit2D(art::PtrVector<recob::Hit> hits,
	       int                        color,
	       evdb::View2D*              view);
    void Draw2DEndPointAndSlope(double        x,
				double        y,
				double        slope,
				int           color,
				evdb::View2D* view);
    void Cluster2D(const art::Event& evt, 
		   evdb::View2D*     view,
		   unsigned int      plane);
    void Prong2D(const art::Event& evt,
		 evdb::View2D*     view,
		 unsigned int      plane);
    void DrawProng2D(std::vector<const recob::Prong*>& prong,
		     evdb::View2D*                     view,
		     unsigned int                      plane,
		     int                               id=-999);
    void Vertex2D(const art::Event& evt,
		  evdb::View2D*     view,
		  unsigned int      plane);
    void Event2D(const art::Event& evt,
		 evdb::View2D*     view,
		 unsigned int      plane);

    void SpacePoint(const recob::Prong* prong,
		    int                 id,
		    evdb::View3D*       view);
    void Prong3D(const art::Event& evt,
		 evdb::View3D*     view);
    void Vertex3D(const art::Event& evt,
		  evdb::View3D*     view);
    void Event3D(const art::Event& evt,
		 evdb::View3D*     view);

    void FillTQHisto(const art::Event&    evt, 
		     unsigned int         plane,
		     unsigned int         wire,
		     TH1F*                histo,
		     std::vector<double>& hstart,
             std::vector<double>& hend,
             std::vector<double>& hitamplitudes,
             std::vector<double>& hpeaktimes);

    void FillQHisto(const art::Event& evt,
		    unsigned int      plane,
		    TH1F*             histo);

    int GetRegionOfInterest(int plane,
			    int& minw,
			    int& maxw,
			    int& mint,
			    int& maxt);
    
    void GetChargeSum(int plane,
		      double& charge,
		      double& convcharge);
    
  private:
    void GetClusterOutlines(art::Ptr<recob::Cluster>& c,
			    std::vector<double>&      tpts,
			    std::vector<double>&      wpts,
			    unsigned int              plane);
    int GetWires(const art::Event&            evt,
		 const std::string&           which,
		 art::PtrVector<recob::Wire>& wires);
    int GetHits(const art::Event&           evt,
		const std::string&          which,
		art::PtrVector<recob::Hit>& hits);
    int GetClusters(const art::Event&               evt,
		    const std::string&              which,
		    art::PtrVector<recob::Cluster>& clust);

    // GetProngs will get Prongs and their derived classes using art::View
    int GetProngs(const art::Event&                 evt,
		  const std::string&                which,
		  std::vector<const recob::Prong*>& prong);

    int GetVertices(const art::Event&              evt,
		    const std::string&             which,
		    art::PtrVector<recob::Vertex>& vertex);

    int GetEvents(const art::Event&              evt,
		  const std::string&             which,
		  art::PtrVector<recob::Event>&  event);

    
        
  private:

    std::vector<unsigned int> fBadChannels; ///< list of bad channels in the detector
    std::vector<int> fWireMin;     ///< lowest wire in interesting region for each plane
    std::vector<int> fWireMax;     ///< highest wire in interesting region for each plane
    std::vector<int> fTimeMin;     ///< lowest time in interesting region for each plane
    std::vector<int> fTimeMax;     ///< highest time in interesting region for each plane
    
    std::vector<double> fRawCharge;     ///< Sum of Raw Charge
    std::vector<double> fConvertedCharge;     ///< Sum of Charge Converted using Birks' formula
    
    
  };
}

#endif
////////////////////////////////////////////////////////////////////////
