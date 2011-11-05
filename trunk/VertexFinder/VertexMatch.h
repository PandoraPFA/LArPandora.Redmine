#ifndef VERTEXMATCH_H
#define VERTEXMATCH_H

#include "TMath.h"
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include "art/Framework/Core/EDProducer.h"

namespace vertex {
   
  class VertexMatch : public art::EDProducer {
    
  public:
    
    explicit VertexMatch(fhicl::ParameterSet const& pset); 
    virtual ~VertexMatch();
    void produce(art::Event& evt);
    
  private:
    std::string fVertexModuleLabel; 
    std::string fHoughModuleLabel;
    double fMaxDistance;
  protected:

    
    
  };


    
}

#endif // VERTEXMATCH_H
