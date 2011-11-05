/// \file    RawDataDrawer.h
/// \brief   Class to aid in the rendering of RawData objects
/// \author  messier@indiana.edu
/// \version $Id: RawDataDrawer.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
#ifndef EVD_RECOBASEDRAWER_H
#define EVD_RECOVASEDRAWER_H
#include <vector>
#ifndef __CINT__
#include "art/Persistency/Common/PtrVector.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#endif

class TH1F;
namespace art  { class Event;     }
namespace evdb { class View2D;    }
namespace evdb { class View3D;    }
namespace raw  { class RawDigit;  }

namespace evd {
  /// Aid in the rendering of RawData objects
  class RawDataDrawer {
  public:
    RawDataDrawer();
    ~RawDataDrawer();

    void RawDigit2D(const art::Event& evt,
		    evdb::View2D*     view,
		    unsigned int      plane);

/*     void RawDigit3D(const art::Event& evt, */
/* 		    evdb::View3D*     view); */

    void FillQHisto(const art::Event& evt, 
		    unsigned int      plane,
		    TH1F*             histo);

    void FillTQHisto(const art::Event& evt, 
		     unsigned int      plane,
		     unsigned int      wire,
		     TH1F*             histo);

    double TotalClockTicks() const { return fTicks; }
    
   
    int GetRegionOfInterest(int plane,
			    int& minw,
			    int& maxw,
			    int& mint,
			    int& maxt);
    
    void GetChargeSum(int plane,
		      double& charge,
		      double& convcharge);
   
  private:
    int GetRawDigits(const art::Event&              evt,
		     art::PtrVector<raw::RawDigit>& rawhits);

    double fTicks;                           ///< number of ticks of the clock
    std::vector<unsigned int> fBadChannels;  ///< bad channels in the detector
   
    std::vector<int> fWireMin;     ///< lowest wire in interesting region for each plane
    std::vector<int> fWireMax;     ///< highest wire in interesting region for each plane
    std::vector<int> fTimeMin;     ///< lowest time in interesting region for each plane
    std::vector<int> fTimeMax;     ///< highest time in interesting region for each plane
    
    std::vector<double> fRawCharge;     ///< Sum of Raw Charge
    std::vector<double> fConvertedCharge;     ///< Sum of Charge Converted using Birks' formula
   };
};

#endif
////////////////////////////////////////////////////////////////////////
