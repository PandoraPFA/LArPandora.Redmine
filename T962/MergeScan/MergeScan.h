////////////////////////////////////////////////////////////////////////
//
// MergeScan Class
// Merge scan information into event record,
// if this info. exists for a given event.
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGESCAN_H
#define MERGESCAN_H

#include "art/Framework/Core/EDProducer.h"
#include "RawData/raw.h"
#include "T962/T962_Objects/ScanInfo.h"

#include <vector>
#include <string>
#include <time.h>

///Merge different data streams 

class TH1F;
class TH2F;

namespace t962 {

  class MergeScan : public art::EDProducer {

  public:
          
    explicit MergeScan(fhicl::ParameterSet const& pset); 
    virtual ~MergeScan();
    void produce(art::Event& evt);
    void beginJob();

  private:

    art::Handle<raw::DAQHeader> fdaq;

  protected: 

    ///<parameters to set
   std::string  daq_modulelabel;               ///< folder for input 
   std::string scanners;
   std::string file;	
   
    bool foundscaninfo;

    	 
  }; // class MergeScan

}

#endif // MERGESCAN_H
