/////////////////////////////////////////////////////////////////////
//
// Definition of track data from MINOSTrackMatch
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <vector>

namespace t962{

   class MINOSTrackMatch  {
      
   public:
      
      MINOSTrackMatch();
      ~MINOSTrackMatch();
      


      int fMINOStrackid;
      int fArgoNeuTtrackid;

     // friend std::ostream& operator << (std::ostream& o, const MINOS& m);

 
   };   
  
}

////////////////////////////////////////////////////////////////////////
