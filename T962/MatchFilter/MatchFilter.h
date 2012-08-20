////////////////////////////////////////////////////////////////////////
//
// MatchFilter class:
// Algoritm to produce a filtered event file having
// events with atleast 1 matched track with MINOS
//
// saima@ksu.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MATCHFILTER_H
#define MATCHFILTER_H

#include "art/Framework/Core/EDFilter.h"
#include "TH1D.h"
#include <string>

namespace filt {

  class MatchFilter : public art::EDFilter  {
    
  public:
    
    explicit MatchFilter(fhicl::ParameterSet const& pset); 
    virtual ~MatchFilter();
         
    
    bool filter(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& p);
    void beginJob();
   

  private: 
 
     std::string fTracks_label;
     std::string fTrackMatchModuleLabel;
     TH1D* fSelectedEvents;

  protected: 

    
  }; // class MatchFilter

}

#endif // MATCHFILTER_H
