#ifndef AGGEVENT_H
#define AGGEVENT_H



#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "AggregateEvent/AggVertex.h"
#include "AggregateEvent/AggTrack.h"
#include "TH2F.h"
#include "TF1.h"
#include <vector>
#include <string>

namespace aggr {
   
  class AggEvent  {
    
  public:
    
    AggEvent(edm::PtrVector<aggr::AggVertex> &,edm::PtrVector<aggr::AggTrack> &);
    ~AggEvent();
    
    
    edm::PtrVector<const recob::Hit> hitlist;
    edm::PtrVector<const recob::Track> track3Dlist;
    edm::PtrVector<const recob::Cluster> clusterlist;
    edm::PtrVector<const recob::Cluster> hclusterlist;
    edm::PtrVector<const recob::Vertex> vertexlist;
    edm::PtrVector<const recob::Shower> showerlist;

    // Next two take C-like structure single argument.
    void SetEventProperties(); // run, event, MC/data, time/date, comprssed, uncprssd size
    void SetCodeProperties(); // ups version, LArSoft version, FPGA version
    // unsigned ints
    void SetTrackCount(unsigned int n)  {fNtracks=n;}
    void SetVertexCount(unsigned int n) {fNvertices=n;}
    void SetShowerCount(unsigned int n) {fNshowers=n;}
    // doubles
    void SetUnCalHitEnergy();
    void SetTotalTrackEnergy();
    void SetTotalShowerEnergy();
    // 4-vector of doubles
    void SetStartPoint_xyzt(); // most upstream hit or vtx?
    void SetEndPoint_xyzt(); // most upstream hit or vtx?

    // Each one of these below classes is created upstream in a 
    // corresponding produce() module (and it has its own
    // pair of *.class files to define it).
      //protected: 
  private:
    edm::PtrVector<aggr::AggVertex> verts; // vertices with associated 3dtracks
    edm::PtrVector<aggr::AggTrack> tracks; // 3dtracks with associated 3dshowers
    unsigned int fNtracks;
    unsigned int fNvertices;
    unsigned int fNshowers;


    /*
    // Don't think I any longer need the friends under the ART port.
    friend class AggVertex;
    friend class AggTrack;
    */

  }; // class AggEvent

}

#endif // AGGEVENT_H
