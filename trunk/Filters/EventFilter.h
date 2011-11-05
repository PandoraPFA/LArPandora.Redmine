////////////////////////////////////////////////////////////////////////
//
// EventFilter class:
//
// This class provides methods for filtering particular events by 
// their evt number.
//  
//
// echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef EVENTFILTER_H
#define EVENTFILTER_H

#include "art/Framework/Core/EDFilter.h"

///filters for events, etc
namespace filter {

  class EventFilter : public art::EDFilter  {

  public:
    
    explicit EventFilter(fhicl::ParameterSet const& ); 
    virtual ~EventFilter();      
    
    bool filter(art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    std::vector < unsigned int > SetOfBadEvents()   const { return fBadEvents;}
    std::vector < unsigned int > SetOfBadRuns()     const { return fBadRuns;  }
    
  private:
    
    std::vector < unsigned int >            fBadEvents; ///< list of bad events
    std::vector < unsigned int >            fBadRuns;   ///< list of bad runs


      }; //class EventFilter
}
#endif // EVENTFILTER_H


