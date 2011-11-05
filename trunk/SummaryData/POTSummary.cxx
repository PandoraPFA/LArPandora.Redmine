////////////////////////////////////////////////////////////////////////
/// \file POTSummary.cxx
/// 
/// Definition of object to store pot related information
/// 
/// \version $Id: RunData.h,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////


#include "SummaryData/POTSummary.h"

namespace sumdata{

  POTSummary::POTSummary():
    totpot(0.0),
    totgoodpot(0.0),
    totspills(0),
    goodspills(0)
  {}
  
  POTSummary::~POTSummary()
  {}
  
  void POTSummary::Print(std::ostream &stream)
  {
    stream<<"This sub run has "<< totspills <<" total spills with an exposure of "
	  <<totpot<<" POT"<<std::endl
	  <<"with cuts on beam quality, there are "<< goodspills 
	  <<" good spills with an exposure of "<<totgoodpot<<std::endl;
  }

}// end namespace
