
#ifndef AGGREGATEVTXANA_H
#define AGGREGATEVTXANA_H

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
#include "FWCore/Framework/interface/EDAnalyzer.h"

// LArSoft includes
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "AggregateEvent/AggVertex.h"

#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include <vector>
#include <string>


namespace aggr {


  class AggregateVertexAna : edm::EDAnalyzer 
  {

  public:

    explicit AggregateVertexAna(edm::ParameterSet const& pset);
    ~AggregateVertexAna();

    void analyze (const edm::Event& evt, edm::EventSetup const&);
    void beginJob(const edm::EventSetup&); 

  private:

    TH1F* HnTrksVtx;
    TH1F* HnVtxes;
    TH1F* HVtxSep;
    TH2F* HVtxRZ;

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fHitModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fVertexModuleLabel;
    std::string fAggVertexModuleLabel;

    edm::PtrVector<recob::Hit> hitlist;
    edm::PtrVector<recob::Vertex> vertexlist;
    edm::PtrVector<recob::Track> tracklist;
    edm::PtrVector<aggr::AggVertex> aggVertexlist;


  }; // class AggregateVertexAna

}  // Namespace aggr

#endif // AGGREGATEVTXANA_H
