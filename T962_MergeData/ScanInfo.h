/////////////////////////////////////////////////////////////////////
// $Id: ScanInfo.h,v 1.3 2010/03/26 19:36:42 brebel Exp $
//  information about the hand scan
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef SCANINFO_H
#define SCANINFO_H

#include <iostream>
#include <fstream>
#include <time.h>
#include <string.h>
#include <iosfwd>

namespace merge {
  
  
  class ScanInfo  {
    
  public:
    ScanInfo();
    ScanInfo(int isneutrino, 
	     int isnotneutrino, 
	     int maybeneutrino, 
	     int trackind, 
	     int trackcol, 
	     int vertindtime, 
	     int vertcoltime, 
	     int vertindwire, 
	     int vertcolwire, 
	     int vnumshower, 
	     int scanner);
    ~ScanInfo();
    
    
    int Get_IsNeutrino() const;
    int Get_IsnotNeutrino() const;
    int Get_IsMaybeNeutrino() const;
    int Get_TrackInd() const;
    int Get_TrackCol() const;
    int Get_VertIndTime() const;
    int Get_VertColTime() const;
    int Get_VertIndWire() const;
    int Get_VertColWire() const;
    int Get_NumShower() const;
    int Get_Scanner() const;

    void SetIsNeutrino(int val);
    void SetIsnotNeutrino(int val);
    void SetIsMaybeNeutrino(int val);
    void SetTrackInd(int val);
    void SetTrackCol(int val);
    void SetVertIndTime(int val);
    void SetVertColTime(int val);
    void SetVertIndWire(int val);
    void SetVertColWire(int val);
    void SetNumShower(int val);
    void SetScanner(int val);

    
  private:
    int isneutrino;
    int isnotneutrino;
    int ismaybeneutrino;
    int trackind;
    int trackcol;
    int vertindtime;
    int vertcoltime;
    int vertindwire;
    int vertcolwire;
    int numshower;
    int scanner;
          
  };
  
  //Non-member functions and operators
  std::ostream& operator<<(std::ostream& , const ScanInfo& );
  
}

////////////////////////////////////////////////////////////////////////
#endif // SCANINFO_H
