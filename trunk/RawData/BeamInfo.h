/////////////////////////////////////////////////////////////////////
// $Id: BeamInfo.h,v 1.3 2010/03/26 19:36:42 brebel Exp $
//  information about the neutrino beam
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef RAWDATA_BEAMINFO_H
#define RAWDATA_BEAMINFO_H

#include <iostream>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iosfwd>

namespace raw {
  
  
  class BeamInfo  {
    
  public:
    BeamInfo();
    BeamInfo(double tor101,double tortgt, double trtgtd,long long int t_ms);
    ~BeamInfo();
    
    
    double get_tor101() const;
    double get_tortgt() const;
    double get_trtgtd() const;
    long long int get_t_ms() const;
    
    
    void SetTOR101(double val);      
    void SetTORTGT(double val);
    void SetTRTGTD(double val);
    void SetT_MS( long long int val);
    
  private:
    double tor101,tortgt,trtgtd;
    long long int t_ms;
    
      
  };
  
  //Non-member functions and operators
  std::ostream& operator<<(std::ostream& , const BeamInfo& );
  
}

////////////////////////////////////////////////////////////////////////
#endif // RAWDATA_BEAMINFO_H
