/////////////////////////////////////////////////////////////////////
// \version $Id: BeamInfo.cxx,v 1.3 2010/03/26 19:36:42 brebel Exp $
// \file    BeamInfo class
// \author  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "RawData/BeamInfo.h"
#include <string.h>
#include <time.h>
#include <cmath>
#include <iomanip>

namespace raw{

  BeamInfo::BeamInfo()
    
  {
    tor101=0;
    tortgt=0;
    trtgtd=0;
    t_ms=0;
    
  }
  
  BeamInfo::~BeamInfo(){ }

  BeamInfo::BeamInfo(double vtor101,double vtortgt, double vtrtgtd,long long int vt_ms)
  {
    tor101=vtor101;
    tortgt=vtortgt;
    trtgtd=vtrtgtd;
    t_ms=vt_ms;
    
  }
  
  
  
  double BeamInfo::get_tor101() const
  {return tor101;}
  
  double BeamInfo::get_tortgt() const
  {return tortgt;}

  double BeamInfo::get_trtgtd() const
  {return trtgtd;}
  
  long long int BeamInfo::get_t_ms() const
  {return t_ms;}

  void  BeamInfo::SetTOR101(double val)
  { tor101=val;}

  void  BeamInfo::SetTORTGT(double val)
  { tortgt=val;}

  void  BeamInfo::SetTRTGTD(double val)
  { trtgtd=val;}

  void  BeamInfo::SetT_MS(long long int val)
  { t_ms=val;}


  std::ostream& operator<<( std::ostream& os, const raw::BeamInfo& o ){//output operator
  
  
    os<<o.get_tor101()<<"\t"<<o.get_tortgt()<<"\t"<<o.get_trtgtd()<<"\t"<<o.get_t_ms();
    
    return os;
  }
}







