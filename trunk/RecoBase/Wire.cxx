////////////////////////////////////////////////////////////////////////
// $Id: Wire.cxx,v 1.10 2010/04/15 18:13:36 brebel Exp $
//
// Wire class
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include "RecoBase/Wire.h"
#include <string>
#include <iostream>
#include <cassert>

namespace recob{

  //----------------------------------------------------------------------
  Wire::Wire() : 
    fSignal(0)
  {

  }

  //----------------------------------------------------------------------
  Wire::Wire(std::vector<double> siglist,
	     art::Ptr<raw::RawDigit> &rawdigit) :
    fSignal(siglist),
    fRawDigit(rawdigit)
  { 

    // uncompress the raw digit adc vector
    std::vector<short> uncompressed(rawdigit->Samples());
    raw::Uncompress(rawdigit->fADC, uncompressed, rawdigit->Compression());

    ///put the pedestal subtracted values of the raw adc's into siglist
    ///if siglist from initializer is empty
    if(fSignal.size()==0)  
      for(unsigned int i = 0; i < uncompressed.size(); ++i)
	fSignal[i] = 1.*uncompressed[i] - rawdigit->GetPedestal();
  
    uncompressed.clear();

  }

  //----------------------------------------------------------------------
  Wire::~Wire()
  {

  }

  //--------------------------------------------------
  double Wire::Signal(int i) const
  {
    char msg[256];
    sprintf(msg,"Wire::ADC(%d) out of range!",i);
    std::cout << msg << std::endl;
    unsigned int j = i;
    assert(i>=0 && (j<fSignal.size()));

    return fSignal[j];
  }

  //--------------------------------------------------
  void Wire::SetSignal(int i, double isig) {
    unsigned int j = i;
    if (fSignal.size()<j+1) fSignal.resize(j+1);
    fSignal[j]=isig;
    return;
  }

}
////////////////////////////////////////////////////////////////////////

