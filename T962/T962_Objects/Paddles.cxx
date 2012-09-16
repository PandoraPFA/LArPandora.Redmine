#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <iomanip>

#include "T962/T962_Objects/Paddles.h"

t962::Paddles::Paddles()
{
  time=0;
  for(int i=0;i<4;++i){
    tdc.resize(4);
    for (int j=0;j<4;++j){
      
      tdc[i].push_back(0);
    }
  }
  
}

t962::Paddles::~Paddles(){}


t962::Paddles::Paddles(time_t t, int pmt1[4], int pmt2[4], int pmt3[4], int pmt4[4])
{
   time=t;
   for (int j=0;j<4;++j)
   {
      tdc[0][j]= pmt1[j];
      tdc[1][j]= pmt2[j];
      tdc[2][j]= pmt3[j];
      tdc[3][j]= pmt4[j];
   }
   
}

time_t t962::Paddles::gettime() const
{return time;}



int t962::Paddles::gettdc(int pmt_number, int hit_number) const
{
  return tdc[pmt_number][hit_number];
}

void t962::Paddles::SetTime(time_t val){
   time=val;
}

void t962::Paddles::SetPMT(int i, int hit[4]){
   for (int j=0;j<4;++j)
      tdc[i][j]= hit[j];
   
}

std::ostream& t962::operator<<( std::ostream& os, const t962::Paddles& o ){//output operator
  
  os<< "time    " << o.gettime();
  std::cout<<std::endl;
  os << std::setprecision(9);
  for (int i=0;i<4;++i){
    os<< "pmt" << i+1 << "   ";
    os<< std::setw(11) << std::right << o.gettdc(i,0)
      << std::setw(11) << std::right << o.gettdc(i,1)
      << std::setw(11) << std::right << o.gettdc(i,2)
      << std::setw(11) << std::right << o.gettdc(i,3);
    std::cout<<std::endl;}
  
  return os;
}




