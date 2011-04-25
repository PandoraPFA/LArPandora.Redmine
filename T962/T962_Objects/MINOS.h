/////////////////////////////////////////////////////////////////////
//
// Definition of track data from MINOS
//
//  kinga.partyka@yale.edu
//  mitchell.soderberg@yale.edu
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <vector>

namespace t962{

   class MINOS  {
      
   public:
      
      MINOS();
      ~MINOS();
      
      void SetMatched(int matched){fmatched = matched;}
      int matched(){return fmatched;}
         
      int frun;
      int fsubRun;
      int fsnarl;
      double futc;
      double fday;
      float ftrkIndex;
      float ftrkChi2;
      double futc1;
      double fnearns;
      float fnearsec;
      float ftrkE;
      float fshwE;
      double fcrateT0;
      double ftmframe;
      double fyear;
      double foffset; 
      float fdtnear;
      float ftrkErange;
      double fsgate53;
      float ftrkqp;
      float ftrkeqp;
      float ftrkVtxX;
      float ftrkVtxY;
      float ftrkVtxZ;
      double ftrkVtxeX;
      double ftrkVtxeY;
      float fcharge;
      float ftrkmom;
      float ftrkVtxT;
      float ftrkTimeT0;
      float ftrkdcosx;
      float ftrkdcosy;
      float ftrkdcosz;
      double fmonth;
      float ftrtgtd;
      float ftortgt;
      float ftor101;
      float ftr101d;
      int fgoodspill;
      int ftrkcontained;
      int fgoodbeam;
      int fntrkstp;
      std::vector<float> ftrkstpX;
      std::vector<float> ftrkstpY;
      std::vector<float> ftrkstpZ;
      std::vector<float> ftrkstpU;
      std::vector<float> ftrkstpV;
               
      //matching
      void SetMINOStrackid(double trackid){fMINOStrackid = trackid;}
      void SetArgoNeuTtrackid(double trackid){fArgoNeuTtrackid = trackid;}
            
      //MC only
      void SetmcPx(double mcpx){fmcPx = mcpx;}
      void SetmcPy(double mcpy){fmcPy = mcpy;}
      void SetmcPz(double mcpz){fmcPz = mcpz;}

      friend std::ostream& operator << (std::ostream& o, const MINOS& m);

   private:    
 
      int fMINOStrackid;
      int fArgoNeuTtrackid; 
      double fmcPx;
      double fmcPy;
      double fmcPz;     
      int fmatched; //our variable to store whether MINOS track is matched to ArgoNeuT track
 
   };   
  
}

////////////////////////////////////////////////////////////////////////
