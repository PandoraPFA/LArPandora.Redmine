////////////////////////////////////////////////////////////////////////
//
// ScanFilter class:
// Tells the downstream rreconstruction to not process events that
// do not pass certain hand scan criteria
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef SCANFILTER_H
#define SCANFILTER_H

#include "art/Framework/Core/EDFilter.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH2D.h"

namespace filt {

  class ScanFilter : public art::EDFilter  {
    
  public:
    
    explicit ScanFilter(fhicl::ParameterSet const& ); 
    virtual ~ScanFilter();
         
    
    bool filter(art::Event& evt);
   
  private: 
 
    std::string fScanModuleLabel;
    int fNeutrino_req, fMinShowers_req, fMinTracks_req, fMaxShowers_req, fMaxTracks_req;
  
  protected: 
    
  }; // class ScanFilter

}

#endif // SCANFILTER_H
