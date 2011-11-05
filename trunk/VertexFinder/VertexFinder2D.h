////////////////////////////////////////////////////////////////////////
/// \file  PrimaryVertexFinder.h
/// \brief Module to find vertices based on 2-d clusters
///
/// \tjyang@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef VertexFinder2D_H
#define VertexFinder2D_H

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include "RecoBase/recobase.h"

class TH1D;

///vertex reconstruction
namespace vertex {
   
 class VertexFinder2D :  public art::EDProducer {
    
  public:
    
    explicit VertexFinder2D(fhicl::ParameterSet const& pset); 
    virtual ~VertexFinder2D();        
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);

  private:
    
    TH1D *dtIC;
  
    std::string fClusterModuleLabel;

  };
    
}



#endif // VertexFinder2D_H
