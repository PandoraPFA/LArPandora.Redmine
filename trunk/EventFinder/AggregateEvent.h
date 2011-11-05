////////////////////////////////////////////////////////////////////////
/// \file  AggregateEvent.h
/// \brief Module to find recob::Event objects
///
/// echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef AGGREGATEEVENT_H
#define AGGREGATEEVENT_H

#include "art/Framework/Core/EDProducer.h"

#include "RecoBase/recobase.h"

#include "TH2F.h"
#include "TF1.h"
#include <vector>
#include <string>

namespace event {
   
  class AggregateEvent : public art::EDProducer {
    
  public:
    
    explicit AggregateEvent(fhicl::ParameterSet const& );
    ~AggregateEvent();

    void produce(art::Event& evt);
    void beginJob();
    
  private:

    std::string fVertexModuleLabel;

    art::PtrVector<recob::Vertex> fvertexlist;

  }; // class AggregateEvent

}

#endif // AGGREGATEEVENT_H
