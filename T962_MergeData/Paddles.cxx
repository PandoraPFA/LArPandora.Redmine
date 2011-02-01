/////////////////////////////////////////////////////////////////////
// \version $Id: Paddles.cxx,v 1.2 2010/02/15 19:34:20 brebel Exp $
// \file    Paddles class
// \author  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "T962_MergeData/Paddles.h"
#include <time.h>
#include <cmath>


merge::Paddles::Paddles()
{
  time=0;
  for(int i=0;i<4;++i){
    tdc.resize(4);
    for (int j=0;j<4;++j){
      
      tdc[i].push_back(0);
    }
  }
  
}

merge::Paddles::~Paddles(){}


merge::Paddles::Paddles(time_t t, int pmt1[4], int pmt2[4], int pmt3[4], int pmt4[4])
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

time_t merge::Paddles::gettime() const
{return time;}



double merge::Paddles::gettdc(int pmt_number, int hit_number) const
{
  return tdc[pmt_number][hit_number];
}

void   merge::Paddles::SetTime(time_t val){
  time=val;
}

void merge::Paddles::SetPMT(int i, int hit[4]){
  for (int j=0;j<4;++j)
    tdc[i][j]= hit[j];
  
}

std::ostream& merge::operator<<( std::ostream& os, const merge::Paddles& o ){//output operator
  
  os<<o.gettime();
  std::cout<<std::endl;
  for (int i=0;i<4;++i){
    
    os<<o.gettdc(i,0)<<'\t'<<o.gettdc(i,1)<<'\t'<<o.gettdc(i,2)<<'\t'<<o.gettdc(i,3);
    std::cout<<std::endl;}
  
  return os;
}




