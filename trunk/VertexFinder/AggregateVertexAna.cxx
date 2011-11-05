////////////////////////////////////////////////////////////////////////
// $Id: AggVertex.cxx,v 1.1 2010/09/02 17:25:11 echurch Exp $
//
// AggregateVertexAna class
//
//  This class will produce the already-discovered vertices
//  containing a vector of track pointers  which we associate
//  to them. 
//
// echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "VertexFinder/AggregateVertexAna.h"
#include "TH1.h"

namespace vertex{

  //-----------------------------------------------
  AggregateVertexAna::AggregateVertexAna(fhicl::ParameterSet const& pset) : 
    fHitModuleLabel(pset.get< std::string >("FFFTHitModuleLabel")),
    fTrack3DModuleLabel(pset.get< std::string >("Track3DModuleLabel")),
    fEndPointModuleLabel(pset.get< std::string >("EndPointModuleLabel")),
    fVertexModuleLabel(pset.get< std::string >("VertexModuleLabel"))
  {


  }

  //-----------------------------------------------
  AggregateVertexAna::~AggregateVertexAna()
  {
  }

  //-----------------------------------------------
  void AggregateVertexAna::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    HnVtxes   = tfs->make<TH1F>("Num Vertices","Num Vertices",8,-0.5,7.5);
    HVtxSep   = tfs->make<TH1F>("Vertices spacing","Vertices spacing",20,0.001,5.0);
    HVtxRZ    = tfs->make<TH2F>("Vtx in RZ","Vtx in RZ",20,-50.0,+50.0,20,0.0,50.0);
    HnTrksVtx = tfs->make<TH1F>("Tracks per vtx","Tracks per vtx",8,-0.5,7.5);

    return;
  }

  //-----------------------------------------------
  void AggregateVertexAna::analyze(const art::Event& evt) 
  {
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    evt.getByLabel(fHitModuleLabel,hitListHandle);
    for(unsigned int ii = 0; ii < hitListHandle->size(); ++ii){
      art::Ptr<recob::Hit> hit(hitListHandle, ii);
      fhitlist.push_back(hit); // class member
    }

    art::Handle< std::vector<recob::EndPoint2D> > epListHandle;
    evt.getByLabel(fEndPointModuleLabel,epListHandle);
    for(unsigned int ii = 0; ii < epListHandle->size(); ++ii){
      art::Ptr<recob::EndPoint2D> ep(epListHandle, ii);
      feplist.push_back(ep); // class member
    }

    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrack3DModuleLabel,trackListHandle);
    for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii){
      art::Ptr<recob::Track> track(trackListHandle, ii);
      ftracklist.push_back(track); // class member
    }

    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii){
      art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      fVertexlist.push_back(vertex); // class member
    }

    HnVtxes->Fill(feplist.size(),1);  

    art::PtrVector<recob::Vertex>::const_iterator avIter = fVertexlist.begin();
    art::PtrVector<recob::Vertex>::const_iterator avIter2 = fVertexlist.begin();
    avIter2++;
    while (avIter != fVertexlist.end())  {            

      art::PtrVector<recob::Track> tvlist = (*avIter)->Tracks();
      HnTrksVtx->Fill(tvlist.size(),1);
      
      if(tvlist.size() < 1) continue;

      art::PtrVector<recob::Hit> bv;
      art::PtrVector<recob::Hit> hitvlist;
      
      // Hits no longer has XYZ() method. To get 3d hit position info I'm going to have to 
      // loop on all the SpacePoints and loop on all Hits from there till it matches
      // one from this vertex. This affects the two Fill() efforts below. EC, 19-Nov-2010.
      
      art::ServiceHandle<geo::Geometry> geom;
      for(unsigned int t = 0;t < geom->NTPC(); ++t){
	for(unsigned int ii = 0; ii < geom->TPC(t).Nplanes(); ++ii){
	  geo::View_t view = geom->Plane(ii).View();
	  bv = tvlist[0]->Clusters(view)[0]->Hits();
	  for (unsigned int jj = 0; jj<bv.size(); jj++) {hitvlist.push_back(bv[jj]);}
	}
      }
      
      art::PtrVector<recob::Hit>::const_iterator hitv = hitvlist.begin();
      //      HVtxRZ->Fill(((TVector3 )(hitv->XYZ())).Perp(),hitv->XYZ()[3],1);
      
      while (avIter2 != fVertexlist.end() && avIter2>avIter)  {            

	art::PtrVector<recob::Hit> bv2;
	art::PtrVector<recob::Hit> hitvlist2;
	for(unsigned int t = 0; t < geom->NTPC(); ++t){
	  for(unsigned int ii = 0; ii < geom->TPC(t).Nplanes(); ++ii){
	    geo::View_t view = geom->Plane(ii).View();
	    bv2 = tvlist[0]->Clusters(view)[0]->Hits();
	    for (size_t jj = 0; jj<bv2.size(); jj++) {hitvlist2.push_back(bv2[jj]);}
	  }
	}
	art::PtrVector<recob::Hit>::const_iterator hitv2 = hitvlist2.begin();
	
	// These two whiles should be each precisely one iteration long.
	while(hitv != hitvlist.end()){
	  while(hitv2 != hitvlist2.end()){
	    //	  TVector3 dist = (TVector3 )(hitv2->XYZ()) - (TVector3 )(hitv->XYZ());
	    TVector3 dist;
	    std::cout << "AggregateVertexAna: dist is " << dist.Mag() << "." << std::endl;
	    HVtxSep->Fill(dist.Mag(),1);
	    hitv2++;
	  }
	  hitv++;
	}
	avIter2++;
      }
      
      avIter++;
    }
    
  }

}// end namespace
