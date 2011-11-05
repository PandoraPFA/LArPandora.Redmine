////////////////////////////////////////////////////////////////////////
// $Id: Track3Dreco.cxx,v 1.36 2011/04/11  Exp $
//
// \file Track3Dreco.h
//
// \author mitch
//
////////////////////////////////////////////////////////////////////////
#ifndef TRACK3DRECO_H
#define TRACK3DRECO_H

#include "art/Framework/Core/EDProducer.h"

#include <vector>
#include <string>

namespace trkf {
   
  class Track3Dreco : public art::EDProducer {
    
  public:
    
    explicit Track3Dreco(fhicl::ParameterSet const& pset);
    ~Track3Dreco();
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:
        
    int             ftmatch; // tolerance for time matching (in time samples) 
    double          fchi2dof;// tolerance for chi2/dof of cluster fit to function
    std::string     fClusterModuleLabel;// label for input cluster collection

  protected: 
    
  
  }; // class Track3Dreco

}

#endif // TRACK3DRECO_H
