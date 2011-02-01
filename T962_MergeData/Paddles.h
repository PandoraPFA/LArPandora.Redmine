/////////////////////////////////////////////////////////////////////
// $Id: Paddles.h,v 1.3 2010/03/26 19:36:42 brebel Exp $
//
// Definition of paddle merge data
//
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef T962_MERGEDATA_PADDLES_H
#define T962_MERGEDATA_PADDLES_H

#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iosfwd>

namespace merge{
  
  class Paddles {
    
  public:
    Paddles();
    Paddles(time_t t, int pmt1[4], int pmt2[4], int pmt3[4], int pmt4[4]);
    ~Paddles();

    time_t gettime() const;
    double gettdc(int pmt_number, int hit_number) const;
    
    void   SetTime(time_t val);      
    void   SetPMT(int i,int hit[4]);
    
  private:
    time_t time;
    std::vector< std::vector<int> > tdc;

  };

  //Non-member functions and operators
  std::ostream& operator<<(std::ostream& , const Paddles& );
  
}

////////////////////////////////////////////////////////////////////////
#endif // MERGEDATA_PADDLES_H
