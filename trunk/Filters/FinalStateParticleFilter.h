////////////////////////////////////////////////////////////////////////
//
// FinalStateParticleFilter class:
// Algoritm to produce a filtered event file having
// events with user-defined final state particles 
//
// saima@ksu.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef FINALSTATEPARTICLEFILTER_H
#define FINALSTATEPARTICLEFILTER_H

#include "art/Framework/Core/EDFilter.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH2D.h"

namespace filt {

  class FinalStateParticleFilter : public art::EDFilter  {
    
  public:
    
    explicit FinalStateParticleFilter(fhicl::ParameterSet const& ); 
    virtual ~FinalStateParticleFilter();
         
    
    bool filter(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& p);
    void beginJob();
   

  private: 
 
    std::string fGenieModuleLabel;
    std::vector<int> fPDG;  
    std::vector<int> fStatusCode;
    TH1D* fSelectedEvents;
    TH1D* fTotalEvents;

  protected: 

    bool isSubset(std::vector<int>& a, std::vector<int>& b);
    
  }; // class FinalStateParticleFilter

}

#endif // FINALSTATEPARTICLEFILTER_H
