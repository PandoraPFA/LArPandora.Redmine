/////////////////////////////////////////////////////////////////////
//
/// \version $Id: SingleGen.h,v 1.2 2010/02/15 19:10:40 brebel Exp $
//  \file  MINOS class
//  \author mitchell.soderberg@yale.edu
////////////////////////////////////////////////////////////////////

#include "T962_MergeData/MINOS.h"

// FOR NEWEST MINOS FILE: ***
merge::MINOS::MINOS(int run, int subRun, int snarl, double utc, double day, float trkIndex,
float trkE, float shwE,float crateT0, float tmframe, double year, std::vector<float> vtx, 
float trkErange,float sgate53, float trkqp, std::vector<float> trkVtx, std::vector<float> trkdcos, 
double month, float trkmom,float charge, float trkstpX[100000], float trkstpY[100000],float trkstpZ[100000],float trkstpU[100000], float trkstpV[100000],float trkeqp,
std::vector<double> trkVtxe,int goodspill,std::vector<float> dtnear_nearsec,std::vector<double> nearns_offset, double utc1,int matched) :
  frun(run), fsubRun(subRun), fsnarl(snarl), futc(utc), fday(day), ftrkIndex(trkIndex),
  ftrkE(trkE),fshwE(shwE), fcrateT0(crateT0), ftmframe(tmframe), fyear(year), fvtx(vtx),
  ftrkErange(trkErange), fsgate53(sgate53), ftrkqp(trkqp), ftrkVtx(trkVtx),
  ftrkdcos(trkdcos),fmonth(month),ftrkmom(trkmom),
  fcharge(charge), ftrkeqp(trkeqp), fgoodspill(goodspill),fdtnear_fnearsec(dtnear_nearsec), fnearns_foffset(nearns_offset),futc1(utc1),fmatched(matched)
{
}
// int goodspill, float dtnear, double nearns, float nearsec, double offset, 
//fgoodspill(goodspill), fdtnear(dtnear), fnearns(nearns), fnearsec(nearsec),foffset(offset),


//fcharge(charge),ftrkstp(trkstp),

// FOR OLD MINOS FILE ***:

// merge::MINOS::MINOS(int run, int subRun, int snarl, double utc, double day, float trkIndex, float trkE, float shwE,	
//                       double crateT0, double tmframe, double year, float vtxX, float vtxY, float vtxZ, float trkErange,	
//                       double sgate53, float trkqp, float trkVtxX, float trkVtxY, float trkVtxZ, float trkdcosx,	
//                       float trkdcosy, float trkdcosz, double month,int matched):	
//   frun(run), fsubRun(subRun), fsnarl(snarl), futc(utc), fday(day), ftrkIndex(trkIndex), ftrkE(trkE),
// 	
//   fshwE(shwE), fcrateT0(crateT0), ftmframe(tmframe), fyear(year), fvtxX(vtxX), fvtxY(vtxY), fvtxZ(vtxZ),
// 	
//   ftrkErange(trkErange), fsgate53(sgate53), ftrkqp(trkqp), ftrkVtxX(trkVtxX), ftrkVtxY(trkVtxY), ftrkVtxZ(trkVtxZ),
// 	
//   ftrkdcosx(trkdcosx), ftrkdcosy(trkdcosy), ftrkdcosz(trkdcosz), fmonth(month),fmatched(matched)
// {
// }

merge::MINOS::MINOS(){}

merge::MINOS::~MINOS(){}



   






std::ostream& merge::operator<<( std::ostream& os, const merge::MINOS& o ){//output operator
  
  std::cout<<std::endl;
  
  return os;

}




