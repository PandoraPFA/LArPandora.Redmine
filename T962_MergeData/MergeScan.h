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

#include "FWCore/Framework/interface/EDProducer.h"
#include "RawData/raw.h"
#include "T962_MergeData/ScanInfo.h"

#include <vector>
#include <string>
#include <time.h>


///Merge different data streams 
class TH1F;
class TH2F;



namespace merge {

  class MergeScan : public edm::EDProducer {

  public:
          
    explicit MergeScan(edm::ParameterSet const& pset); 
    virtual ~MergeScan();
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);

  private:

  protected: 

    ///<parameters to set
   std::string  daq_modulelabel;               ///< folder for input 
   std::vector<std::string> scanners;
   std::string file;	
   
    bool foundscaninfo;

    	 
  }; // class MergeScan

}

#endif // MERGESCAN_H
