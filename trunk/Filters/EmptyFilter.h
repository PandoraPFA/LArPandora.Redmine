////////////////////////////////////////////////////////////////////////
//
// EmptyFilter class:
// Algorith to produce event files with the
// blank events removed using only hit information.
//
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef EMPTYFILTER_H
#define EMPTYFILTER_H

#include "art/Framework/Core/EDFilter.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH2D.h"

namespace filt {

  class EmptyFilter : public art::EDFilter  {
    
  public:
    
    explicit EmptyFilter(fhicl::ParameterSet const& ); 
    virtual ~EmptyFilter();
         
    
    bool filter(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& p);
    void beginJob();
   

  private: 
 
    std::string fHitsModuleLabel;
    double  fMinIonization;  
    int fMinNumHits;
    TH1I * totHitHist;
    TH1I * selHitHist;
    TH1I * rejHitHist;
    TH2D * totIonSelHist;
    TH2D * totIonRejHist;
    TH1I * numEventHist;
    TH2I * resultTable;  

  protected: 
    
  }; // class EmptyFilter

}

#endif // EMPTYFILTER_H
