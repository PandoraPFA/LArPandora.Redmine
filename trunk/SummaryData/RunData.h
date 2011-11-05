////////////////////////////////////////////////////////////////////////
/// \file RunData.h
/// 
/// Definition of object to store run related information
/// 
/// \version $Id: RunData.h,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SD_RUNDATA_H
#define SD_RUNDATA_H
#include "Geometry/geo.h"

namespace sumdata {

  class RunData
  {

  public: 
    
    RunData(); // Default constructor
    explicit RunData(geo::DetId_t detid);
    ~RunData();

    geo::DetId_t DetId() const { return fDetId; }

  private:

    geo::DetId_t fDetId;  ///< detector id 

  };
}

#endif // SD_RUNDATA_H
