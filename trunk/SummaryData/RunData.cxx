////////////////////////////////////////////////////////////////////////
/// \file RunData.cxx
/// 
/// Definition of object to store run related information
/// 
/// \version $Id: RunData.cxx,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SummaryData/RunData.h"

namespace sumdata {

  //---------------------------------------------------------
  RunData::RunData() 
    : fDetId(geo::kArgoNeuT)
  {
  }

  //---------------------------------------------------------
  RunData::RunData(geo::DetId_t detid) 
    : fDetId(detid)
  {
  }

  //---------------------------------------------------------
  RunData::~RunData()
  {
  }

}  
