/////////////////////////////////////////////////////////////////////
//
// Definition of paddle data
//
//  kinga.partyka@yale.edu
//  msoderbe@syr.edu
////////////////////////////////////////////////////////////////////
#ifndef PADDLES_H
#define PADDLES_H

#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iosfwd>

namespace t962{
  
  class Paddles {
    
  public:
    Paddles();
    Paddles(time_t t, int pmt1[4], int pmt2[4], int pmt3[4], int pmt4[4]);
    ~Paddles();

    time_t gettime() const;
    int    gettdc(int pmt_number, int hit_number) const;
    
    void   SetTime(time_t val);      
    void   SetPMT(int i,int hit[4]);
    
  private:
    time_t time;
    std::vector< std::vector<int> > tdc;

  };

  //Non-member functions and operators
  std::ostream& operator<<(std::ostream& , const Paddles& );
  
}//namespace t962

////////////////////////////////////////////////////////////////////////
#endif // PADDLES_H
