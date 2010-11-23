#ifndef AGGREGATEEVENT_H
#define AGGREGATEEVENT_H

#include "FWCore/Framework/interface/EDProducer.h"

#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "AggregateEvent/AggVertex.h"
#include "TH2F.h"
#include "TF1.h"
#include <vector>
#include <string>

namespace aggr {
   
  class AggregateEvent : public edm::EDProducer {
    
  public:
    
    explicit AggregateEvent(edm::ParameterSet const& );
    ~AggregateEvent();

    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob();
    
  private:

    std::string fHitModuleLabel;
    std::string fT3DModuleLabel;
    std::string fDBScanModuleLabel;
    std::string fHCModuleLabel;
    std::string fVertexModuleLabel;
    std::string fShowerModuleLabel;
    std::string fCalorimetryModuleLabel;


    edm::PtrVector<recob::Cluster> clusterlist;
    edm::PtrVector<recob::Cluster> hclusterlist;
    edm::PtrVector<recob::Hit> hitlist;
    edm::PtrVector<recob::Vertex> vertexlist;
    edm::PtrVector<recob::Track> tracklist;
    edm::PtrVector<recob::Shower> showerlist;
    unsigned int runNum;
    unsigned int evtNum;
  }; // class AggregateEvent

}

#endif // AGGREGATEEVENT_H
