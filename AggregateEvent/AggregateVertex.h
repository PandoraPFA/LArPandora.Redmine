
#ifndef AGGREGATEVTX_H
#define AGGREGATEVTX_H

// Framework includes
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/View.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
#include "FWCore/Framework/interface/EDProducer.h" 

// LArSoft includes
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include <vector>
#include <string>

namespace edm {
  class Event;
  class ParameterSet;
}
 
namespace aggr {


  class AggregateVertex : public edm::EDProducer
  {

  public:

    explicit AggregateVertex(edm::ParameterSet const& pset);
    virtual ~AggregateVertex();

    void produce(edm::Event& evt, edm::EventSetup const&); 
    void beginJob(const edm::EventSetup&); 

    //std::vector<aggr::AggVertex *> MatchV2T(edm::PtrVector<recob::Vertex> , edm::PtrVector<recob::Track>);    
    std::auto_ptr< std::vector<aggr::AggVertex> >  MatchV2T();

  private:

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fVertexModuleLabel;

    edm::PtrVector<recob::Cluster> clusterlist;
    edm::PtrVector<recob::Cluster> hclusterlist;
    edm::PtrVector<recob::Hit> hitlist;
    edm::PtrVector<recob::Vertex> vertexlist;
    edm::PtrVector<recob::Track> tracklist;
    edm::PtrVector<recob::Shower> showerlist;
    edm::PtrVector<recob::Vertex> vertexlistStrong;

    const recob::Vertex* copiedVert;
    std::vector<const recob::Track*> matchedTracks;

  }; // class AggregateVertex

}  // Namespace aggr

#endif // AGGREGATEVTX_H
