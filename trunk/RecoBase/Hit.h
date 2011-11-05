////////////////////////////////////////////////////////////////////////////
// \version $Id: Hit.h,v 1.7 2010/03/01 21:32:44 bpage Exp $
//
// \brief Definition of Hit object for LArSoft
//
// \author mitchell.soderberg@yale.edu
//
////////////////////////////////////////////////////////////////////////////
#ifndef HIT_H
#define HIT_H

#include <vector>
#include <iosfwd>
#include <string>
#include <iostream>
#include <iomanip>

#include "art/Persistency/Common/Ptr.h"

#include "RecoBase/Wire.h"
#include "Geometry/geo.h"

namespace recob {

  ///hits are 2D representations of charge deposited in the tdc/wire plane
  ///hits are assumed to be made from deconvoluted, unipolar signals
  class Hit {
    public:
      Hit(); // Default constructor
      explicit Hit(art::Ptr<recob::Wire> &wire, 
		   double startTime, double sigmaStartTime,
		   double endTime,   double sigmaEndTime, 
		   double peakTime,  double sigmaPeakTime,
		   double totcharge, double sigmaTotCharge, 
		   double maxcharge, double sigmaMaxCharge, 
		   int    multiplicity,
		   double goodnessOfFit);
      
      ~Hit();

      // Get Methods
      double              StartTime()                 const { return fStartTime;     }	       	 
      double              EndTime()        	      const { return fEndTime;       }	       	 
      double              PeakTime()       	      const { return fPeakTime;      }	       	 
      double              SigmaStartTime() 	      const { return fSigmaStartTime;}	       	 
      double              SigmaEndTime()   	      const { return fSigmaEndTime;  }	       	 
      double              SigmaPeakTime()  	      const { return fSigmaPeakTime; }	       	 
      int                 Multiplicity()   	      const { return fMultiplicity;  }	       	 
      unsigned short      Channel()        	      const { return Wire()->RawDigit()->Channel();}
      double              Charge(bool max=false)      const;
      double              SigmaCharge(bool max=false) const;
      double              GoodnessOfFit()    	      const { return fGoodnessOfFit; }
      
      geo::View_t           View()                    const { return fView;}
      art::Ptr<recob::Wire> Wire()                    const { return fWire;}

      std::vector<double> fHitSignal;//vector of ADC values within the hit window
      
      friend std::ostream& operator << (std::ostream & o, const Hit & a);
      friend bool          operator <  (const Hit & a, const Hit & b);

private:
      
      double                fStartTime;      ///< initial tdc tick for hit
      double                fSigmaStartTime; ///< uncertainty on initial tick
      double                fEndTime;        ///< final tdc tick for hit
      double                fSigmaEndTime;   ///< uncertainty on final tick
      double                fPeakTime;       ///< tdc for the peak charge deposition
      double                fSigmaPeakTime;  ///< uncertainty for tdc of the peak
      double                fCharge;         ///< total charge deposited for hit
      double                fMaxCharge;      ///< maximum ADC value in hit window
      double                fSigmaCharge;    ///< uncertainty in total charge deposited
      double                fSigmaMaxCharge; ///< maximum ADC value in hit window
      int                   fMultiplicity;   ///< how many hits could this one be shared with
      double                fGoodnessOfFit;  ///< how well do we believe we know this hit?
      art::Ptr<recob::Wire> fWire;           ///< index of Wire object this Hit was created on
      geo::View_t           fView;           ///< view for hit
    
    };
}

#endif
