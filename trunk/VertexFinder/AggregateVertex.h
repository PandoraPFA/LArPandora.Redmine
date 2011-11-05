////////////////////////////////////////////////////////////////////////
/// \file  AggregateVertex.h
/// \brief Module to find vertices based on 2-d clusters
///
/// echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef AGGREGATEVERTEX_H
#define AGGREGATEVERTEX_H

// Framework includes
#include "art/Framework/Core/EDProducer.h" 

// LArSoft includes
#include "RecoBase/recobase.h"

#include <vector>
#include <string>

namespace vertex {


  class AggregateVertex : public art::EDProducer
  {

  public:

    explicit AggregateVertex(fhicl::ParameterSet const& pset);
    virtual ~AggregateVertex();

    void produce(art::Event& evt); 
    void beginJob(); 

    std::auto_ptr< std::vector<recob::Vertex> >  MatchV2T();

  private:

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fEndPointModuleLabel;

    art::PtrVector<recob::EndPoint2D> feplist;
    art::PtrVector<recob::Track>      ftracklist;
    art::PtrVector<recob::EndPoint2D> feplistStrong;

  }; // class AggregateVertex

}  // Namespace vertex

#endif // AGGREGATEVERTEX_H
