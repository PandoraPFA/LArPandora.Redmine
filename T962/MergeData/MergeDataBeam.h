////////////////////////////////////////////////////////////////////////
//
// Merge POT information into DAQ480 data,
// if this info. exists for a given event.
// Kinga.partyka@yale.edu
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGEDATABEAM_H
#define MERGEDATABEAM_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Persistency/Common/Ptr.h"

#include "RawData/BeamInfo.h"
#include "RawData/DAQHeader.h"

#include <string>

///T962 analysis for merging different types of data
namespace merge {

  class MergeDataBeam : public art::EDProducer {

  public:
          
    explicit MergeDataBeam(fhicl::ParameterSet const& pset); 
    virtual ~MergeDataBeam();

    void produce(art::Event& evt);
  
  private: 

    bool MergeBeam(raw::BeamInfo& beam);                 ///method to merge beam data
    
    art::Handle<raw::DAQHeader> fdaq;
    std::string  fdaq_modulelabel;               ///< folder for input 
     
    
  }; // class MergeDataBeam

}

#endif // MERGEDATABEAM_H
