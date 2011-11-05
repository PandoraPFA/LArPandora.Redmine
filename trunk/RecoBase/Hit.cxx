////////////////////////////////////////////////////////////////////////
// \version $Id: Hit.cxx,v 1.7 2010/02/15 20:32:46 brebel Exp $
//
// \brief Definition of Hit reconstruction object
//
// \author mitchell.soderberg@yale.edu
////////////////////////////////////////////////////////////////////////

#include "RecoBase/Hit.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace recob{

  //----------------------------------------------------------------------
  Hit::Hit() 
    : fStartTime(0.)
    , fSigmaStartTime(0.)
    , fEndTime(0.)
    , fSigmaEndTime(0.)
    , fPeakTime(0.) 
    , fSigmaPeakTime(0.)
    , fCharge(0.)
    , fMaxCharge(0.)
    , fSigmaCharge(0.)
    , fSigmaMaxCharge(0.)
    , fMultiplicity(1)
    , fGoodnessOfFit(0.)
    , fView(geo::kUnknown)
  {
    mf::LogWarning("RecoBase") << "you are using the default recob::Hit() construction - "
			       << "that's a bad idea" << std::endl;
    fHitSignal.clear();
  }

  //----------------------------------------------------------------------
  Hit::Hit(art::Ptr<recob::Wire> &wire,
	   double startTime, double sigmaStartTime,
	   double endTime,   double sigmaEndTime, 
	   double peakTime,  double sigmaPeakTime,
	   double totcharge, double sigmaTotCharge, 
	   double maxcharge, double sigmaMaxCharge, 
	   int    multiplicity,
	   double goodnessOfFit) :
    fStartTime(startTime),
    fSigmaStartTime(sigmaStartTime),
    fEndTime(endTime),
    fSigmaEndTime(sigmaEndTime),
    fPeakTime(peakTime), 
    fSigmaPeakTime(sigmaPeakTime),
    fCharge(totcharge),
    fMaxCharge(maxcharge),
    fSigmaCharge(sigmaTotCharge),
    fSigmaMaxCharge(sigmaMaxCharge),
    fMultiplicity(multiplicity),
    fGoodnessOfFit(goodnessOfFit),
    fWire(wire)
  {

    art::ServiceHandle<geo::Geometry> geo;
    unsigned int p = 0; 
    unsigned int w = 0;
    unsigned int t = 0;
    geo->ChannelToWire(wire->RawDigit()->Channel(), t, p, w);

    fView = geo->Plane(p,t).View();

    fHitSignal.clear();
  }

  //----------------------------------------------------------------------
  Hit::~Hit()
  {
  }

  //----------------------------------------------------------------------
  double Hit::Charge(bool max) const
  {
    if(max) return fMaxCharge;
    return fCharge;
  }

  //----------------------------------------------------------------------
  double Hit::SigmaCharge(bool max) const
  {
    if(max) return fSigmaMaxCharge;
    return fSigmaCharge;
  }

  //----------------------------------------------------------------------
  // ostream operator.  
  //

  std::ostream& operator<< (std::ostream & o, const Hit & a)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << " Channel "        << std::setw(5) << std::right << a.Channel() 
      << " View = "         << std::setw(3) << std::right << a.View() 
      << " StartTime = "    << std::setw(7) << std::right << a.StartTime() 
      << " +/- "            << std::setw(7) << std::right << a.SigmaStartTime() 
      << " EndTime = "      << std::setw(7) << std::right << a.EndTime() 
      << " +/- "            << std::setw(7) << std::right << a.SigmaEndTime() 
      << " PeakTime = "     << std::setw(7) << std::right << a.PeakTime()  
      << " Charge = "       << std::setw(7) << std::right << a.Charge()
      << " +/- "            << std::setw(7) << std::right << a.SigmaCharge() 
      << " Multiplicity = " << std::setw(5) << std::right << a.Multiplicity() 
      << " GoodnessOfFit = "<< std::setw(7) << std::right << a.GoodnessOfFit() ;

    return o;
  }


  //----------------------------------------------------------------------
  // < operator.  
  //
  bool operator < (const Hit & a, const Hit & b)
  {
    if(a.Channel() != b. Channel())
      return a.Channel()<b.Channel();
    if(a.View() != b.View())
      return a.View()<b.View();
    if(a.StartTime() != b.StartTime())
      return a.StartTime() < b.StartTime();

    return false; //They are equal
  }


  //----------------------------------------------------------------------
}
