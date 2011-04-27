/////////////////////////////////////////////////////////////////////
//
//  MINOS class 
//
////////////////////////////////////////////////////////////////////

#include "T962/T962_Objects/MINOS.h"

namespace t962{

   MINOS::MINOS():
      fmatched(-1),
      ftrkstpX(1,-999.0),
      ftrkstpY(1,-999.0),
      ftrkstpZ(1,-999.0),
      ftrkstpU(1,-999.0),
      ftrkstpV(1,-999.0)
   {
   }
   
   MINOS::~MINOS()
   {
   }

   //--------------------------------------------------------------------
   // ostream operator.
   //
   std::ostream& operator << (std::ostream& o, const MINOS& m)
   {
      o << resetiosflags(std::ios::fixed);
      o << " MINOS TrkIndex "  << std::setw(2) << std::right << m.ftrkIndex 
        << "  : Q = " << std::setw(3) << std::right << m.fcharge;
      o << std::setiosflags(std::ios::fixed) << std::setprecision(3);
      o << " TrkVtx = (" << std::setw(7) << std::right << m.ftrkVtxX
        << " " << std::setw(7) << std::right << m.ftrkVtxY
        << " " << std::setw(7) << std::right << m.ftrkVtxZ << " )"
        << " StartCosines = (" << std::setw(7) << std::right << m.ftrkdcosx
        << " " << std::setw(7) << std::right << m.ftrkdcosy
        << " " << std::setw(7) << std::right << m.ftrkdcosz << " )"
        << std::setprecision(2) 
        << " E = " << std::setw(5) << std::right << m.ftrkE
        << " Erange = " << std::setw(5) << std::right << m.ftrkErange
        << " Mom = " << std::setw(7) << std::right << m.ftrkmom
        << " Chi2 = " << std::setw(7) << std::right << m.ftrkChi2;

      return o;
   }

}






