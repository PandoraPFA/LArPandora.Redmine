////////////////////////////////////////////////////////////////////////
// $Id: RawDigit.cxx,v 1.13 2010/03/26 19:36:42 brebel Exp $
//
// RawDigit class
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include "RawData/RawDigit.h"
#include <string>
#include <iostream>
#include <cassert>

namespace raw{

  //----------------------------------------------------------------------
  RawDigit::RawDigit() : 
    fADC(0), 
    fChannel(0), 
    fSamples(0), 
    fPedestal(0.), 
    fSigma(0.),
    fCompression(raw::kNone)
  {

  }

  //----------------------------------------------------------------------
  RawDigit::RawDigit(unsigned short channel,
		     unsigned short samples,
		     std::vector<short> adclist, 
		     raw::Compress_t compression) :
    fADC(adclist), 
    fChannel(channel), 
    fSamples(samples),
    fPedestal(0.), 
    fSigma(0.),
    fCompression(compression)
  { 

  }

  //----------------------------------------------------------------------
  RawDigit::RawDigit(unsigned int channel,
		     std::vector<short> adclist,
		     raw::Compress_t compression) :
    fADC(adclist), 
    fChannel(channel), 
    fSamples(0),
    fPedestal(0.), 
    fSigma(0.),
    fCompression(compression)
  {

  }

  //----------------------------------------------------------------------
  RawDigit::~RawDigit()
  {

  }

  //--------------------------------------------------
  short RawDigit::ADC(int i) const
  {
    unsigned int j = i;
    assert(i>=0 && (j<fADC.size()));

    return fADC[j];
  }

  //--------------------------------------------------
  void RawDigit::SetADC(int i, short iADC) 
  {
    unsigned int j = i;
    if (fADC.size()<j+1) fADC.resize(j+1);
    fADC[j]=iADC;
    return;
  }

  //----------------------------------------------------------------------
  void RawDigit::Set(unsigned int channel, int i, short adc)
  {
    fChannel = channel;
    unsigned int j = i;
    if (fADC.size()<j+1) fADC.resize(j+1);
    fADC[j] = adc;
  }

  //----------------------------------------------------------------------
  void RawDigit::SetPedestal(double ped)
  {

    fPedestal = ped;
    fSigma = 1.;

  }
}
////////////////////////////////////////////////////////////////////////

