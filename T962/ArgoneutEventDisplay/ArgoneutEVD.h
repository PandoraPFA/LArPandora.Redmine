////////////////////////////////////////////////////////////////////////
//
// The module to control the event display
//
// brebel@fnal.gov
// msoderbe@syr.edu
////////////////////////////////////////////////////////////////////////
#ifndef ARGONEUTEVD_H
#define ARGONEUTEVD_H
#ifndef __CINT__

// Framework Includes
#include "art/Framework/Core/EDAnalyzer.h"

#include <string>
#include "TH1D.h"

/// The Event Display
namespace argoevd{

  /// a class for transporting photons in a roughly realistic way
  class ArgoneutEVD : public art::EDAnalyzer 
  {
   public:
     explicit ArgoneutEVD(fhicl::ParameterSet const &pset);
     virtual ~ArgoneutEVD();

     void analyze(art::Event const& evt);
     void beginJob();

   };
}
#endif // __CINT__
#endif // EVD
////////////////////////////////////////////////////////////////////////
