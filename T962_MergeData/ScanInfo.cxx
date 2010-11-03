/////////////////////////////////////////////////////////////////////
// \version $Id: ScanInfo.cxx,v 1.3 2010/03/26 19:36:42 brebel Exp $
// \file    ScanInfo class
// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "T962_MergeData/ScanInfo.h"
#include <string.h>
#include <time.h>
#include <cmath>
#include <iomanip>

namespace merge{

  ScanInfo::ScanInfo()
    
  {
    isneutrino=-1;
    isnotneutrino=-1;
    ismaybeneutrino=-1;
    trackind=-1;
    trackcol=-1;
    vertindtime=-1;
    vertcoltime=-1;
    vertindwire=-1;
    vertcolwire=-1;
    numshower=-1;
  }
  
  ScanInfo::~ScanInfo(){ }

  ScanInfo::ScanInfo(int visneutrino,int visnotneutrino,int vismaybeneutrino, int vtrackind, int vtrackcol, int vvertindtime,int vvertcoltime,int vvertindwire,int vvertcolwire, int vnumshower)
  {
    isneutrino=visneutrino;
    isnotneutrino=visnotneutrino;
    ismaybeneutrino=vismaybeneutrino;
    trackind=vtrackind;
    trackcol=vtrackcol;
    vertindtime=vvertindtime;
    vertcoltime=vvertcoltime;
    vertindwire=vvertindwire;
    vertcolwire=vvertcolwire;
    numshower=vnumshower;
  }
  
  
  
  int ScanInfo::Get_IsNeutrino() const
  {return isneutrino;}
  
  int ScanInfo::Get_IsnotNeutrino() const
  {return isnotneutrino;}
  
  int ScanInfo::Get_IsMaybeNeutrino() const
  {return ismaybeneutrino;}
  
  int ScanInfo::Get_TrackInd() const
  {return trackind;}
  
  int ScanInfo::Get_TrackCol() const
  {return trackcol;}
  
  int ScanInfo::Get_VertIndTime() const
  {return vertindtime;}
    
  int ScanInfo::Get_VertColTime() const
  {return vertcoltime;}
  
  int ScanInfo::Get_VertIndWire() const
  {return vertindwire;}
    
  int ScanInfo::Get_VertColWire() const
  {return vertcolwire;}
  
  int ScanInfo::Get_NumShower() const
  {return numshower;}


  void ScanInfo::SetIsNeutrino(int val) 
  {isneutrino=val;}
  
  void ScanInfo::SetIsnotNeutrino(int val) 
  {isnotneutrino=val;}

  void ScanInfo::SetIsMaybeNeutrino(int val) 
  {ismaybeneutrino=val;}
    void ScanInfo::SetTrackInd(int val) 
  {trackind=val;}
  
  void ScanInfo::SetTrackCol(int val) 
  {trackcol=val;}
  
  void ScanInfo::SetVertIndTime(int val) 
  {vertindtime=val;}
  
  void ScanInfo::SetVertColTime(int val) 
  {vertcoltime=val;}  
  
  void ScanInfo::SetVertIndWire(int val) 
  {vertindwire=val;}
  
  void ScanInfo::SetVertColWire(int val) 
  {vertcolwire=val;}  

   void ScanInfo::SetNumShower(int val) 
  {numshower=val;} 


  std::ostream& operator<<( std::ostream& os, const merge::ScanInfo& o ){//output operator
  
  
    os<<o.Get_IsnotNeutrino()<<"\t"<<o.Get_IsMaybeNeutrino()<<"\t"<<o.Get_IsNeutrino()<<"\t"<<o.Get_TrackInd()<<"\t"<<o.Get_TrackCol()<<"\t"<<o.Get_VertIndTime()<<"\t"<<o.Get_VertColTime()<<o.Get_VertIndWire()<<"\t"<<o.Get_VertColWire()<<"\t"<<o.Get_NumShower();
    
    return os;
  }
}







