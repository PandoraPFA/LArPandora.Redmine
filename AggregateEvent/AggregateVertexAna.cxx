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

#include "AggregateEvent/AggregateVertexAna.h"
#include "AggregateEvent/AggVertex.h"
#include "TH1.h"

//-----------------------------------------------
aggr::AggregateVertexAna::AggregateVertexAna(edm::ParameterSet const& pset) : 
  fHitModuleLabel(pset.getParameter< std::string >("FFFTHitModuleLabel")),
  fTrack3DModuleLabel(pset.getParameter< std::string >("Track3DModuleLabel")),
  fVertexModuleLabel(pset.getParameter< std::string >("VertexModuleLabel")),
  fAggVertexModuleLabel(pset.getParameter< std::string >("AggVertexModuleLabel"))
{


  HnVtxes   = new TH1F("Num Vertices","Num Vertices",8,-0.5,7.5);
  // Below is always one. Don't make it.
  //HnhinVtxes = new TH1F("Num hits in Vertices","Num Hits in Vertices",5,-0.5,4.5);
  HVtxSep = new TH1F("Vertices spacing","Vertices spacing",20,0.001,5.0);
  HVtxRZ = new TH2F("Vtx in RZ","Vtx in RZ",20,-50.0,+50.0,20,0.0,50.0);
  HnTrksVtx = new TH1F("Tracks per vtx","Tracks per vtx",8,-0.5,7.5);

}

void aggr::AggregateVertexAna::analyze(const edm::Event& evt, edm::EventSetup const&) 
{
  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);
  for(unsigned int ii = 0; ii < hitListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Hit> hit(hitListHandle, ii);
      hitlist.push_back(hit); // class member
    }
  edm::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      vertexlist.push_back(vertex); // class member
    }
  edm::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fTrack3DModuleLabel,trackListHandle);
  for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Track> track(trackListHandle, ii);
      tracklist.push_back(track); // class member
    }
  edm::Handle< std::vector<aggr::AggVertex> > aggVertexListHandle;
  evt.getByLabel(fAggVertexModuleLabel,aggVertexListHandle);
  for(unsigned int ii = 0; ii < aggVertexListHandle->size(); ++ii)
    {
      edm::Ptr<aggr::AggVertex> aggVertex(aggVertexListHandle, ii);
      aggVertexlist.push_back(aggVertex); // class member
    }

  HnVtxes->Fill(vertexlist.size(),1);  


  edm::PtrVectorItr<aggr::AggVertex> avIter = aggVertexlist.begin();
  edm::PtrVectorItr<aggr::AggVertex> avIter2 = aggVertexlist.begin();
  avIter2++;
  while (avIter != aggVertexlist.end())  
    {            

      edm::PtrVector<recob::Track> tvlist = (*avIter)->matchedTracks;
      HnTrksVtx->Fill(tvlist.size(),1);

      edm::Ptr<recob::Vertex> vtx = (*avIter)->copiedVert;
      edm::PtrVector<recob::Hit> bv;
      edm::PtrVector<recob::Hit> hitvlist;

      // Hits no longer has XYZ() method. To get 3d hit position info I'm going to have to 
      // loop on all the SpacePoints and loop on all Hits from there till it matches
      // one from this vertex. This affects the two Fill() efforts below. EC, 19-Nov-2010.

      edm::Service<geo::Geometry> geom;
      int nplv = geom->Nplanes();
      for(int ii = 0; ii < nplv; ++ii)
	{
	  geo::View_t view = geom->Plane(ii).View();
	  bv = vtx->Hits(view);
	  for (int jj = 0; jj<bv.size(); jj++) {hitvlist.push_back(bv[jj]);}
	}
      edm::PtrVectorItr<recob::Hit> hitv = hitvlist.begin();
      //      HVtxRZ->Fill(((TVector3 )(hitv->XYZ())).Perp(),hitv->XYZ()[3],1);

      while (avIter2 != aggVertexlist.end() && avIter2>avIter)  
	{            
	  edm::Ptr<recob::Vertex> vtx2 = (*avIter2)->copiedVert;
	  edm::PtrVector<recob::Hit> bv2;
	  edm::PtrVector<recob::Hit> hitvlist2;
	  int nplv2 = geom->Nplanes();
	  for(int ii = 0; ii < nplv2; ++ii)
	    {
	      geo::View_t view = geom->Plane(ii).View();
	      bv2 = vtx2->Hits(view);
	  for (int jj = 0; jj<bv2.size(); jj++) {hitvlist2.push_back(bv2[jj]);}
	    }
	  edm::PtrVectorItr<recob::Hit> hitv2 = hitvlist2.begin();

	  // These two whiles should be each precisely one iteration long.
 	  while(hitv != hitvlist.end())
	    {
	      while(hitv2 != hitvlist2.end())
		{
		  //	  TVector3 dist = (TVector3 )(hitv2->XYZ()) - (TVector3 )(hitv->XYZ());
		  TVector3 dist;
		  std::cout << "aggr::AggregateVertexAna: dist is " << dist.Mag() << "." << std::endl;
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
