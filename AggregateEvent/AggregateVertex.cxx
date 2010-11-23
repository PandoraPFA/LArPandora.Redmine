////////////////////////////////////////////////////////////////////////
// $Id: AggregateVertex.cxx,v 1.1 2010/09/02 17:25:11 echurch Exp $
//
// AggregateVertex class
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

#include "AggregateEvent/AggVertex.h"
#include "AggregateEvent/AggregateVertex.h"
#include "TH1.h"

//-----------------------------------------------
aggr::AggregateVertex::AggregateVertex(edm::ParameterSet const& pset) : 
  fDBScanModuleLabel(pset.getParameter< std::string >("DBScanModuleLabel")),
  fHoughModuleLabel(pset.getParameter< std::string >("HoughModuleLabel")),
  fTrack3DModuleLabel(pset.getParameter< std::string >("Track3DModuleLabel")),
  fVertexModuleLabel(pset.getParameter< std::string >("VertexModuleLabel"))
{
  produces< std::vector<aggr::AggVertex> >();
}

void aggr::AggregateVertex::produce(edm::Event& evt, edm::EventSetup const&) 
{

  
  //  std::auto_ptr<std::vector<aggr::AggVertex> > vcol(new std::vector<aggr::AggVertex>);
  /*
    Deviating from usual ART script here in the way in which we put stuff on the event.
    (Namely, I'm not doing it by means of usual above declaration,)
    This is because I want to let a function create the pointer to the new objects,
    AggVertexes, and hand them back to this produce() method. Check this site:
    http://www.gotw.ca/publications/using_auto_ptr_effectively.htm. EC, 21-Nov-2010.
   */


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


  // Only match strong vertices to tracks.
  edm::PtrVectorItr<recob::Vertex> vIter = vertexlist.begin();

  while (vIter != vertexlist.end())  
    {            
      edm::Ptr <recob::Vertex> vtx = (*vIter);  
      if (vtx->ID() < 3) vertexlistStrong.push_back(vtx); // -- cuz there's one less
      vIter++;
    }

  // We will suck up the hits out of each vertex and out of the tracks
  // and see if there's an overlap of a sufficient (>0) number of hits. If so,
  // call the track a match to the vtx. Then, .... stick the track pointer(s)
  // into the AggVertex object.  EC, 23-July-2010.

  std::auto_ptr< std::vector<aggr::AggVertex> > vcol (MatchV2T());

  evt.put(vcol);
}

//std::vector<aggr::AggVertex *> aggr::AggregateVertex::MatchV2T()
std::auto_ptr< std::vector<aggr::AggVertex> >  aggr::AggregateVertex::MatchV2T()
{
  std::cout << "aggr::AggregateEvent::MatchV2T(): (strong) vertexlistStrong and tracklist lengths are " << vertexlistStrong.size() << " and " <<  tracklist.size() << "." << std::endl;


  // Bail if there are no tracks or vertices.
  //  if (!((int)vertexlistStrong.size()) || !((int)tracklist.size())) return NULL;
  if (vertexlistStrong.isNull() || tracklist.isNull()) 
    {
      std::vector<aggr::AggVertex *> np; 
      return std::auto_ptr< std::vector<aggr::AggVertex> > (new std::vector<aggr::AggVertex>);
    }
  edm::Service<geo::Geometry> geom;
  // Loop on the vertices, and all the hits in each
  edm::PtrVectorItr<recob::Vertex> vIter = vertexlistStrong.begin();
  //  std::vector<aggr::AggVertex *> verts;
  std::auto_ptr< std::vector<aggr::AggVertex> > verts(new std::vector<aggr::AggVertex>);

  while (vIter != vertexlistStrong.end())  
    {            
      edm::PtrVector<recob::Track> tlistAssoc; // Will fill with matching tracks.
      edm::Ptr <recob::Vertex> vtx = (*vIter);
      edm::PtrVector<recob::Hit> hitvertexlistStrong;
      
      edm::PtrVector<recob::Hit> bv;
      int nplv = geom->Nplanes();
      for(int ii = 0; ii < nplv; ++ii)
	{
	  geo::View_t view = geom->Plane(ii).View();
	  // std::cout << "aggr::AggregateEvent::MatchV2T(): Vtx hit view is " << view << std::endl;
	  bv = vtx->Hits(view);
	  for (int jj=0; jj < bv.size(); jj++) hitvertexlistStrong.push_back(bv[jj]);
	}

      // Should be just one hit per vtx, as per Josh, but we loop anyway.
      edm::PtrVectorItr<recob::Hit> hvIter = hitvertexlistStrong.begin();
      while (hvIter != hitvertexlistStrong.end())  
	{
	  edm::Ptr<recob::Hit> hitv = (*hvIter);
	  // Now loop on the tracks and all the clusters in each, and all
	  // the hits in each of those.
	  edm::PtrVectorItr<recob::Track> tIter = tracklist.begin();
	  
	  while (tIter != tracklist.end())  
	    {
	      
	      edm::Ptr<recob::Track> trk = (*tIter);
	      
	      // Don't use this one, cuz these are Maddolena's
	      // added hits. I want the originals: the U,V hits.
	      // std::vector<const recob::Cluster*>clustlist3 = trk->Clusters(geo::k3D);
	      edm::PtrVector<recob::Cluster>clustlistU = trk->Clusters(geo::kU);
	      edm::PtrVector<recob::Cluster>clustlistV = trk->Clusters(geo::kV);
	      edm::PtrVector<recob::Cluster>clustlist;
	      for (int jj=0; jj < clustlistU.size(); jj++) {clustlist.push_back(clustlistU[jj]);}
	      for (int jj=0; jj < clustlistV.size(); jj++) {clustlist.push_back(clustlistV[jj]);}
	      
	      edm::PtrVectorItr<recob::Cluster> cIter = clustlist.begin();
	      
	      while (cIter != clustlist.end())
		{
		  edm::Ptr<recob::Cluster> cls = (*cIter);
		  
		  edm::PtrVector<recob::Hit> hittlist;
		  edm::PtrVector<recob::Hit> b;
		  

		  int nplc = geom->Nplanes();
		  for(int ii = 0; ii < nplc; ++ii)
		    {
		      geo::View_t view = geom->Plane(ii).View();
		      b = cls->Hits(view);
		      for (int jj=0; jj < b.size(); jj++) hittlist.push_back(b[jj]);
		    }
		  
		  /*
		    In Maddalena's code, they were all k3D hits.
		    First cluster was always 2 hits: the starting/ending 
		    hits. Second cluster was full of hits on 3d track. 
		    But now -- in her code, trkf::Track3Dreco.cxx --
		    I've packed original hits too, which I'll use. So, now
		    First two clusters are orignal hits in original views.
		    Third cluster is the 2 hits: the starting/ending hits.
		    Fourth cluster contains MP's matched hits on 3d track.
		  */

		  if (!hittlist.size()) hittlist = cls->Hits(geo::k3D);
		  
		  edm::PtrVectorItr<recob::Hit> htIter = hittlist.begin();
		  while (htIter != hittlist.end())  
		    {

		      edm::Ptr<recob::Hit> hitt = (*htIter);
		      
		      if (hitt==hitv)
			{
			  //std::cout << "aggr::AggregateEvent::MatchV2T(): BooYiggity! Vtx/Trk hit match! hitt, hitv= " << hitt << " " <<hitv << std::endl;
			  //std::cout << "aggr::AggregateEvent::MatchV2T(): vtx, trk " << vtx<< " " <<trk << std::endl;
			  tlistAssoc.push_back(trk);
			  htIter = hittlist.end()-1; // jump to end of track hitlist, since we've satisfied match
			  cIter = clustlist.end()-1; // jump to end of track cluster list too.
			}			
		      htIter++;
		    } // end hits in clusters in track loop
		  cIter++;
		} // end cluster in track loop
	      tIter++;
	    } // end track loop
	  hvIter++;
	}

      // Now if matching tracks were found for this vertex then create the AggVert object.
      
      if (tlistAssoc.size()>0)
	    {
	      // std::vector<const recob::Track*> tclistAssoc = tlistAssoc;
	      aggr::AggVertex*  vertTmp = new aggr::AggVertex(vtx);
	      vertTmp->SetMatchedTracks(tlistAssoc);
	      verts->push_back(*vertTmp);
	    }
	      
      vIter++;
	  //HnhinVtxes->Fill(hitvertexlistStrong.size(),1);
	      
    }     // end vIter vtx loop 
    

  return verts;

      
} // end AggregateVertex::MatchV2T


