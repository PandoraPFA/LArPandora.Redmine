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
   int frun;
    int fsubRun;
    int fsnarl;
    double futc;
    double fday;
    float ftrkIndex;
    float ftrkE;
    float fshwE;
    float fcrateT0;
    float ftmframe;
    double fyear;
    //float fvtxX;
    //float fvtxY;
    //float fvtxZ;
    float ftrkErange;
    float fsgate53;
    float ftrkqp;
    //float ftrkVtxX;
    //float ftrkVtxY;
    //float ftrkVtxZ;
    //float ftrkdcosx;
    //float ftrkdcosy;
    //float ftrkdcosz;
    double fmonth;
    int fmatched;
    float ftrkmom;
    float fcharge;
    float ftrkstpU[100000];
    float ftrkstpV[100000];
    int fntrkstp;
    float ftrkstpX[100000];
    float ftrkstpY[100000];
    float ftrkstpZ[100000];
    float ftrkeqp;
    // ftrkVtxeX;
    //float ftrkVtxeY;
     std::vector<float> fvtx;
    std::vector<float> ftrkVtx;
    std::vector<float>ftrkdcos;
     //std::vector<float> ftrkstp;
     std::vector<double> ftrkVtxe;
     
     std::vector<float> fdtnear_fnearsec;
     int fgoodspill;
     std::vector<double> fnearns_foffset;
     //float fnearsec;
     //double foffset; 
     double futc1;
     
    
    //*** FOR OLD MINOS FILES:***
    
    // int frun;
//      int fsubRun;
// 	
//     int fsnarl;
// 	
//     double futc;
// 	
//     double fday;
// 	
//     float ftrkIndex;
// 	
//     float ftrkE;
// 	
//     float fshwE;
// 	
//     double fcrateT0;
// 	//
//     double ftmframe;
// 	//
//     double fyear;
// 	
//     float fvtxX;
// 	
//     float fvtxY;
// 	
//     float fvtxZ;
// 	
//     float ftrkErange;
// 	//
//     double fsgate53;
// 	//
//     float ftrkqp;
// 	
//     float ftrkVtxX;
// 	
//     float ftrkVtxY;
// 	
//     float ftrkVtxZ;
// 	
//     float ftrkdcosx;
// 	
//     float ftrkdcosy;
// 	
//     float ftrkdcosz;
// 	
//     double fmonth;
// 	
//     int fmatched;
//     
    
    
  public:
  // FOR NEWEST MINOS FILE: ***
  
      MINOS(int run, int subRun, int snarl, double utc, double day, float trkIndex, 
     float trkE, float shwE,float crateT0, float tmframe, double year, std::vector<float> fvtx,
     float trkErange,float sgate53, float trkqp, std::vector<float> trkVtx, 
     std::vector<float> trkdcos, double month,float trkmom, 
     float charge, float trkstpX[100000], float trkstpY[100000],float trkstpZ[100000],float trkstpU[100000],float trkstpV[100000], float trkeqp,
     std::vector<double> ftrkVtxe,int goodspill,std::vector<float> dtnear_fnearsec,std::vector<double> nearns_foffset,double utc1, int matched);
     
     
     
    // int goodspill,float dtnear,double nearns, float nearsec, double offset, 
     
   //FOR OLD MINOS FILE: ***
   
   // MINOS(int run, int subRun, int snarl, double utc, double day, float trkIndex, float trkE, float shwE,
// 	
//           double crateT0, double tmframe, double year, float vtxX, float vtxY, float vtxZ, float trkErange,
// 	
//           double sgate53, float trkqp, float trkVtxX, float trkVtxY, float trkVtxZ, float trkdcosx,
// 	
//           float trkdcosy, float trkdcosz, double month, int matched);
   
   
   
   MINOS();
    ~MINOS();
    
    //old file:
    //  void SetMatched(int matched)
//     {fmatched=matched;}
// 	void SetRun(int run)
//     {frun=run;}
// 	void SetSubRun(int subrun)
//     {fsubRun=subrun;}
//     void SetSnarl(int snarl)
//     {fsnarl=snarl;}
// 	void SetUtc(double utc)
//     {futc=utc;}
//     void SetDay(double day)
//     {fday=day;}
//     void SetTrkIndex(float trkIndex)
//     {ftrkIndex=trkIndex;}
//     void SetTrkE(float trkE)
//     {ftrkE=trkE;}
//     void SetShwE(float shwE)
//     {fshwE=shwE;}
//     void SetCrateT0(double crateT0)
//     {fcrateT0=crateT0;}
//     void SetTmframe(double tmframe)
//     {ftmframe=tmframe;}
//     void SetYear(double year)
//     {fyear=year;}
//     void SetVtxX(float vtxX)
//     {fvtxX=vtxX;}
//     void SetVtxY(float vtxY)
//     {fvtxY=vtxY;}
//     void SetVtxZ(float vtxZ)
//     {fvtxZ=vtxZ;}
//     void SetTrkErange(float trkErange)
//     {ftrkErange=trkErange;}
//     void SetSgate53(double sgate53)
//     {fsgate53=sgate53;}
//     void SetTrkqp(float trkqp)
//     {ftrkqp=trkqp;}
//     void SetTrkVtxX(float TrkVtxX)
//     {ftrkVtxX=TrkVtxX;}
//      void SetTrkVtxY(float TrkVtxY)
//     {ftrkVtxY=TrkVtxY;}
//      void SetTrkVtxZ(float TrkVtxZ)
//     {ftrkVtxZ=TrkVtxZ;}
//      void SetTrkdcosx(float trkdcosx)
//     {ftrkdcosx=trkdcosx;}
//     void SetTrkdcosy(float trkdcosy)
//     {ftrkdcosy=trkdcosy;}
//     void SetTrkdcosz(float trkdcosz)
//     {ftrkdcosz=trkdcosz;}
//     void SetMonth(double month)
//     {fmonth=month;}
    
     //new file:
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
    void SetVtx(float vtxX,float vtxY,float vtxZ)
    {fvtx[0]=vtxX;
    fvtx[1]=vtxY;
    fvtx[2]=vtxZ;
    }
    
    void SetTrkErange(float trkErange)
    {ftrkErange=trkErange;}
    void SetSgate53(double sgate53)
    {fsgate53=sgate53;}
    void SetTrkqp(float trkqp)
    {ftrkqp=trkqp;}
    void SetTrkVtx(float TrkVtxX,float TrkVtxY,float TrkVtxZ)
    {ftrkVtx[0]=TrkVtxX;
    ftrkVtx[1]=TrkVtxY;
    ftrkVtx[2]=TrkVtxZ;}
    
     void SetTrkdcos(float trkdcosx,float trkdcosy,float trkdcosz)
    {ftrkdcos[0]=trkdcosx;
    ftrkdcos[1]=trkdcosy;
    ftrkdcos[2]=trkdcosz;
    }
    
    void SetTrkmom(float trkmom)
    {ftrkmom=trkmom;
    }
    void SetCharge(float charge)
    {fcharge=charge;
    }
    void SetTrkStpX(float trkstpX[100000])
    {ftrkstpX[100000]=trkstpX[100000];}
    void SetTrkStpY(float trkstpY[100000])
    {ftrkstpY[100000]=trkstpY[100000];}
    void SetTrkStpZ(float trkstpZ[100000])
    {ftrkstpZ[100000]=trkstpZ[100000];}
    void SetTrkStpU(float trkstpU[100000])
    {ftrkstpU[100000]=trkstpU[100000];}
    void SetTrkStpV(float trkstpV[100000])
    {ftrkstpV[100000]=trkstpV[100000];}
    
    
    
    void SetTrkeqp(float trkeqp)
    {ftrkeqp=trkeqp;}
    void SetTrkVtxe(double trkVtxeX, double trkVtxeY)
    {
    ftrkVtxe[0]=trkVtxeX;
    ftrkVtxe[1]=trkVtxeY;
    }
    void SetGoodspill( int goodspill)
    {fgoodspill=goodspill;}
    void SetDtnear_nearsec(float dtnear,float nearsec)
    {fdtnear_fnearsec[0]=dtnear;
    fdtnear_fnearsec[1]=nearsec;}
    void SetNearns_offset(double nearns,double offset)
    {fnearns_foffset[0]=nearns;
    fnearns_foffset[1]=offset;}
    //void SetNearsec(float nearsec)
    //{fnearsec=nearsec;}
   // void SetOffset(double offset)
    //{foffset=offset;}
    void SetUtc1(double utc1)
    {futc1=utc1;}
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
