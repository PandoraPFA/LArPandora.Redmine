////////////////////////////////////////////////////////////////////////
// $Id: SpacePts.cxx,v 1.36 2011/04/11  Exp $
//
// \file SpacePts.h
//
// \author mitch
//
////////////////////////////////////////////////////////////////////////
#ifndef SPACEPTS_H
#define SPACEPTS_H

#include "art/Framework/Core/EDProducer.h" //

#include <vector>
#include <string>

namespace trkf {
   
  class SpacePts : public art::EDProducer {
    
  public:
    
    explicit SpacePts(fhicl::ParameterSet const& pset);
    ~SpacePts();
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:
        
    int             ftmatch; // tolerance for time matching (in time samples) 
    double          fPreSamplings; // in ticks
    double fvertexclusterWindow;
    std::string     fClusterModuleLabel;// label for input cluster collection
    std::string     fEndPoint2DModuleLabel;//label for input EndPoint2D collection
  protected: 
    
  
  }; // class SpacePts

}

#endif // SPACEPTS_H
