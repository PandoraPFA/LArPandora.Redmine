////////////////////////////////////////////////////////////////////////
//
// Merge Paddles information into ArgoNeuT data,
// if this info. exists for a given event.
// 
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGEDATAPADDLES_H
#define MERGEDATAPADDLES_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"

#include "T962/T962_Objects/Paddles.h"
#include "RawData/DAQHeader.h"

#include <string>

///T962 analysis for merging different types of data
namespace merge {

  class MergeDataPaddles : public art::EDProducer {

  public:
          
    explicit MergeDataPaddles(fhicl::ParameterSet const& pset); 
    virtual ~MergeDataPaddles();

    void produce(art::Event& evt);
  
  private: 

    bool MergePaddles(t962::Paddles& paddles, art::Event& evt);   ///method to merge paddles data
    
    art::Handle<raw::DAQHeader> fdaq;
    std::string  fdaq_modulelabel;               ///< folder for input 
     
    
  }; // class MergeDataPaddles

}

#endif // MERGEDATAPADDLES_H
