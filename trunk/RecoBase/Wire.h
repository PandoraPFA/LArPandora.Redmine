////////////////////////////////////////////////////////////////////////
// \version $Id: Wire.h,v 1.7 2010/03/26 20:07:08 brebel Exp $
//
// \brief Definition of basic wire object.  The deconvoluted signals are stored in this class
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef WIRE_H
#define WIRE_H

#include <vector>
#include <iosfwd>
#include "RawData/raw.h"

#include "art/Persistency/Common/Ptr.h"


///Reconstruction base classes
namespace recob {
  
  class Wire {
    public:
      Wire(); // Default constructor
      Wire(std::vector<double> siglist,
	   art::Ptr<raw::RawDigit> &rawdigit);
      
      ~Wire();

      // Set Methods
      void            SetSignal(int i, double isig);  
      void	      SetStatus(unsigned int s)     {fStatus = s;}

      // Get Methods
      double                   Signal(int i) const;      
      unsigned int             NSignal()     const { return fSignal.size(); } 
      unsigned int             Status()      const { return fStatus;}
      art::Ptr<raw::RawDigit>  RawDigit()    const { return fRawDigit;}

      std::vector<double> fSignal;

private:

      unsigned int            fStatus; 
      art::Ptr<raw::RawDigit> fRawDigit; ///vector to index of raw digit for this wire

    };
}

#endif // WIRE_H

////////////////////////////////////////////////////////////////////////
