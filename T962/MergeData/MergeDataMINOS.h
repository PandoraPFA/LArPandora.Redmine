////////////////////////////////////////////////////////////////////////
//
// Merge MINOS information into DAQ480 data,
// if this info. exists for a given event.
// joshua.spitz@yale.edu
// kinga.partyka@yale.edu
// soderber@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGEDATAMINOS_H
#define MERGEDATAMINOS_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Persistency/Common/Ptr.h"

#include "RawData/DAQHeader.h"
#include "RawData/BeamInfo.h"
#include "T962/T962_Objects/MINOS.h"

#include <string>
#include <vector>


class TH1D;
class TH2D;

namespace merge {
   
   class MergeDataMINOS : public art::EDProducer {
      
   public:
      
      explicit MergeDataMINOS(fhicl::ParameterSet const& pset); 
      virtual ~MergeDataMINOS();
      
      void beginJob();
      void produce(art::Event& evt);
      
   private:

      TH2D* fproblemevent2d;
      TH1D* fPOTdiff_matched;
      TH2D* fMINOSrun_event;
      TH1D* futc1_tms_diff;

      
   protected: 
      
      bool MergeMINOS(std::vector<t962::MINOS> & vec_minos); //method to merge MINOS data

      void parse_filename(std::string filename, int & minos_t1, int & minos_t2); //method to grab relevant MINOS info. from filenames

      art::Ptr<raw::DAQHeader> fdaq;
      std::string  fdaq_modulelabel;               // label for input DAQHeader object 
      
      art::Ptr<raw::BeamInfo> fbeam;
      std::string fbeam_modulelabel;               // label for input BeamInfo object
      
   }; // class MergeDataMINOS
   
}

#endif // MERGEDATAMINOS_H
