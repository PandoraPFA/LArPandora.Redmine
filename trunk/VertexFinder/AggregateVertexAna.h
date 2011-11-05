///////////////////////////////////////////////////////////////////////
/// \file    AggregateVertexAna.h
/// \brief   
/// \author  echurch@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////

#ifndef VERTEX_AGGREGATEVERTEXANA_H
#define VERTEX_AGGREGATEVERTEXANA_H

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"

// LArSoft includes
#include "RecoBase/recobase.h"

#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include <vector>
#include <string>


namespace vertex {

  class AggregateVertexAna : public art::EDAnalyzer {

  public:

    explicit AggregateVertexAna(fhicl::ParameterSet const& pset);
    ~AggregateVertexAna();

    void analyze (const art::Event& evt);
    void beginJob(); 

  private:

    TH1F* HnTrksVtx;
    TH1F* HnVtxes;
    TH1F* HVtxSep;
    TH2F* HVtxRZ;

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fHitModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fEndPointModuleLabel;
    std::string fVertexModuleLabel;

    art::PtrVector<recob::Hit>        fhitlist;
    art::PtrVector<recob::EndPoint2D> feplist;
    art::PtrVector<recob::Track>      ftracklist;
    art::PtrVector<recob::Vertex>     fVertexlist;


  }; // class AggregateVertexAna

}  // Namespace aggr

#endif // AGGREGATEVTXANA_H
