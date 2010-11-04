/////////////////////////////////////////////////////////////////////
//
// Definition of track data from MINOS
//
//  kinga.partyka@yale.edu
//  mitchell.soderberg@yale.edu
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iosfwd>
#include <vector>

namespace merge{
  
  
  class MINOS  {

  private:
  //*** for NEWEST MINOS FILE:***
   //  int frun;
//     int fsubRun;
//     int fsnarl;
//     double futc;
//     double fday;
//     float ftrkIndex;
//     float ftrkE;
//     float fshwE;
//     float fcrateT0;
//     float ftmframe;
//     double fyear;
//     //float fvtxX;
//     //float fvtxY;
//     //float fvtxZ;
//     float ftrkErange;
//     float fsgate53;
//     float ftrkqp;
//     //float ftrkVtxX;
//     //float ftrkVtxY;
//     //float ftrkVtxZ;
//     //float ftrkdcosx;
//     //float ftrkdcosy;
//     //float ftrkdcosz;
//     double fmonth;
//     int fmatched;
//     float ftrkmom;
//     int fcharge;
//     // float ftrkstpU;
// //     float ftrkstpV;
//    // int fntrkstp;
//    //  float ftrkstpX;
// //     float ftrkstpY;
// //     float ftrkstpZ;
//     float ftrkeqp;
//     // ftrkVtxeX;
//     //float ftrkVtxeY;
//      std::vector<float> fvtx;
//     std::vector<float> ftrkVtx;
//     std::vector<float>ftrkdcos;
//      std::vector<float> ftrkstp;
//      std::vector<float> ftrkVtxe;
    
    //*** FOR OLD MINOS FILES:***
    
    int frun;
     int fsubRun;
	
    int fsnarl;
	
    double futc;
	
    double fday;
	
    float ftrkIndex;
	
    float ftrkE;
	
    float fshwE;
	
    double fcrateT0;
	//
    double ftmframe;
	//
    double fyear;
	
    float fvtxX;
	
    float fvtxY;
	
    float fvtxZ;
	
    float ftrkErange;
	//
    double fsgate53;
	//
    float ftrkqp;
	
    float ftrkVtxX;
	
    float ftrkVtxY;
	
    float ftrkVtxZ;
	
    float ftrkdcosx;
	
    float ftrkdcosy;
	
    float ftrkdcosz;
	
    double fmonth;
	
    int fmatched;
    
    
    
  public:
  // FOR NEWEST MINOS FILE: ***
  
     //  MINOS(int run, int subRun, int snarl, double utc, double day, float trkIndex, 
//      float trkE, float shwE,float crateT0, float tmframe, double year, std::vector<float> fvtx,
//      float trkErange,float sgate53, float trkqp, std::vector<float> trkVtx, 
//      std::vector<float> trkdcos, double month,float trkmom, 
//      int charge, std::vector<float> trkstp, float trkeqp,
//      std::vector<float> ftrkVtxe, int matched);
     
   //FOR OLD MINOS FILE: ***
   
   MINOS(int run, int subRun, int snarl, double utc, double day, float trkIndex, float trkE, float shwE,
	
          double crateT0, double tmframe, double year, float vtxX, float vtxY, float vtxZ, float trkErange,
	
          double sgate53, float trkqp, float trkVtxX, float trkVtxY, float trkVtxZ, float trkdcosx,
	
          float trkdcosy, float trkdcosz, double month, int matched);
   
   
   
   MINOS();
    ~MINOS();
    
     void SetMatched(int matched)
    {fmatched=matched;}
	
	 void SetRun(int run)
    {frun=run;}
	void SetSubRun(int subrun)
    {fsubRun=subrun;}
    void SetSnarl(int snarl)
    {fsnarl=snarl;}
	void SetUtc(double utc)
    {futc=utc;}
    void SetDay(double day)
    {fday=day;}
    void SetTrkIndex(float trkIndex)
    {ftrkIndex=trkIndex;}
    void SetTrkE(float trkE)
    {ftrkE=trkE;}
    void SetShwE(float shwE)
    {fshwE=shwE;}
    void SetCrateT0(double crateT0)
    {fcrateT0=crateT0;}
    void SetTmframe(double tmframe)
    {ftmframe=tmframe;}
    void SetYear(double year)
    {fyear=year;}
    void SetVtxX(float vtxX)
    {fvtxX=vtxX;}
    void SetVtxY(float vtxY)
    {fvtxY=vtxY;}
    void SetVtxZ(float vtxZ)
    {fvtxZ=vtxZ;}
    void SetTrkErange(float trkErange)
    {ftrkErange=trkErange;}
    void SetSgate53(double sgate53)
    {fsgate53=sgate53;}
    void SetTrkqp(float trkqp)
    {ftrkqp=trkqp;}
    void SetTrkVtxX(float TrkVtxX)
    {ftrkVtxX=TrkVtxX;}
     void SetTrkVtxY(float TrkVtxY)
    {ftrkVtxY=TrkVtxY;}
     void SetTrkVtxZ(float TrkVtxZ)
    {ftrkVtxZ=TrkVtxZ;}
     void SetTrkdcosx(float trkdcosx)
    {ftrkdcosx=trkdcosx;}
    void SetTrkdcosy(float trkdcosy)
    {ftrkdcosy=trkdcosy;}
    void SetTrkdcosz(float trkdcosz)
    {ftrkdcosz=trkdcosz;}
    void SetMonth(double month)
    {fmonth=month;}
    
    
    
    
 int getMatched(){

   return fmatched;

    }
    
  };

  //Non-member functions and operators
  std::ostream& operator<<(std::ostream& , const MINOS& );
  
}

////////////////////////////////////////////////////////////////////////
