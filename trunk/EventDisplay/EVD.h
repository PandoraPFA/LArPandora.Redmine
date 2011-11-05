////////////////////////////////////////////////////////////////////////
// $Id: EVD.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
//
// The module to control the event display
//
// brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVD_EVD_H
#define EVD_EVD_H
#ifndef __CINT__

// Framework Includes
#include "art/Framework/Core/EDAnalyzer.h"

#include <string>
#include "TH1D.h"

/// The Event Display
namespace evd{

  /// a class for transporting photons in a roughly realistic way
  class EVD : public art::EDAnalyzer 
  {
   public:
     explicit EVD(fhicl::ParameterSet const &pset);
     virtual ~EVD();

     void analyze(art::Event const& evt);
     void beginJob();

   };
}
#endif // __CINT__
#endif // EVD
////////////////////////////////////////////////////////////////////////
