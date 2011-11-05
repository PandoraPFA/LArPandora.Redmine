////////////////////////////////////////////////////////////////////////
// $Id: DAQHeader.cxx,v 1.7 2010/02/15 19:34:20 brebel Exp $
//
// DAQHeader class
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include "RawData/DAQHeader.h"
#include <string>
#include <iostream>

namespace raw{

  //----------------------------------------------------------------------
  DAQHeader::DAQHeader() : 
    fStatus(0), fFixed(0), fFormat(0), fSoftware(0), fRun(0), fEvent(0), fTime(0), fSpare(0), fNchan(0)
  {

  }

  //----------------------------------------------------------------------
  DAQHeader::DAQHeader(unsigned int status) : 
    fStatus(status), fFixed(0), fFormat(0), fSoftware(0), fRun(0), fEvent(0), fTime(0), fSpare(0), fNchan(0)
  {

  }

  //----------------------------------------------------------------------
  DAQHeader::DAQHeader(unsigned int status,
		       int fixed,
		       unsigned short format,
		       unsigned short software,
		       unsigned short run,
		       unsigned short event,
		       time_t time,
		       short spare,
		       unsigned short nchan) :
    fStatus(status), fFixed(fixed), fFormat(format), fSoftware(software), fRun(run), fEvent(event), fTime(time), fSpare(spare), fNchan(nchan)
  { 
  
  }

  //----------------------------------------------------------------------
  DAQHeader::~DAQHeader()
  {

  }
}
////////////////////////////////////////////////////////////////////////

