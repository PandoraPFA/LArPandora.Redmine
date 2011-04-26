////////////////////////////////////////////////////////////////////////
//
// Merge MINOS simuation information onto event record
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGESIMMINOS_H
#define MERGESIMMINOS_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Persistency/Common/Ptr.h"

#include "RawData/DAQHeader.h"
#include "RawData/BeamInfo.h"
#include "T962/T962_Objects/MINOS.h"
#include "Simulation/sim.h"

#include <string>
#include <vector>


class TH1D;
class TH2D;

namespace merge {
   
   class MergeSimMINOS : public art::EDProducer {
      
   public:
      
      explicit MergeSimMINOS(fhicl::ParameterSet const& pset); 
      virtual ~MergeSimMINOS();
      
      void beginJob();
      void produce(art::Event& evt);
      
   private:
      
   protected: 
      
      bool MergeMINOS(std::vector<t962::MINOS> & vec_minos); //method to merge MINOS data
      art::Handle< std::vector<sim::Particle> > parHandle;
      art::Handle<raw::BeamInfo> fbeam;
      std::string fbeam_modulelabel;               // label for input BeamInfo object
      std::string  fG4ModuleLabel;
   }; // class MergeSimMINOS
   
}

#endif // MERGESIMMINOS_H
