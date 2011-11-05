////////////////////////////////////////////////////////////////////////
/// \file POTSummary.h
/// 
/// Definition of object to store pot related information
/// 
/// \version $Id: RunData.h,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef POTSUM_H
#define POTSUM_H

#include <iostream>

namespace sumdata {

  class POTSummary {
  public: 
    POTSummary();
    ~POTSummary();

    void Print(std::ostream &stream);

    double totpot;
    double totgoodpot;

    int totspills;
    int goodspills;
  };

}
#endif //POTSUM_H
     

    
    
    

