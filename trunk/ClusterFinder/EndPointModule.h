////////////////////////////////////////////////////////////////////////
/// \file  EndPointModule.h
/// \brief Module to find 2D end points
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef EndPointModule_H
#define EndPointModule_H

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>

///vertex reconstruction
namespace cluster {

  ///module to find 2D end points 
 class EndPointModule :  public art::EDProducer {
    
  public:
    
    explicit EndPointModule(fhicl::ParameterSet const& pset); 
    virtual ~EndPointModule();        

    void produce(art::Event& evt);
    
  private:

    std::string fDBScanModuleLabel;
   
  };
    
}



#endif // EndPointModule_H
