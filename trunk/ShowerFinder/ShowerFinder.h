////////////////////////////////////////////////////////////////////////
/// \file  ShowerFinder.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  roxanne
////////////////////////////////////////////////////////////////////////
#ifndef SHOWERFINDER_H
#define SHOWERFINDER_H

#include <vector>
#include <string>

#include "art/Framework/Core/EDProducer.h" 

///shower finding
namespace shwf {
   
  class ShowerFinder : public art::EDProducer {
    
  public:
    
    explicit ShowerFinder(fhicl::ParameterSet const&); 
    virtual ~ShowerFinder();
         
    void produce(art::Event& evt); 
 

  private:

    std::string     fVertexModuleLabel;  
    std::string     fClusterModuleLabel;  
    std::string     fHoughLineModuleLabel;  
    std::string     fVertexStrengthModuleLabel;  
    double     fRcone;  
    double     fLcone;  

  protected: 
    
    
  }; // class ShowerFinder
  
  
}

#endif // SHOWERFINDER_H
