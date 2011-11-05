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

#include "VertexFinder/AggregateVertex.h"
#include "TH1.h"

namespace vertex {

  //-----------------------------------------------
  AggregateVertex::AggregateVertex(fhicl::ParameterSet const& pset) : 
    fDBScanModuleLabel(pset.get< std::string >("DBScanModuleLabel")),
    fHoughModuleLabel(pset.get< std::string >("HoughModuleLabel")),
    fTrack3DModuleLabel(pset.get< std::string >("Track3DModuleLabel")),
    fEndPointModuleLabel(pset.get< std::string >("EndPointModuleLabel"))
  {
    produces< std::vector<recob::Vertex> >();
  }
  
  //-----------------------------------------------
  AggregateVertex::~AggregateVertex()
  {
  }

  //-----------------------------------------------
  void AggregateVertex::beginJob()
  {
  }

  //-----------------------------------------------
  void AggregateVertex::produce(art::Event& evt) 
  {

  
    //  std::auto_ptr<std::vector<AggVertex> > vcol(new std::vector<AggVertex>);
    /*
      Deviating from usual ART script here in the way in which we put stuff on the event.
      (Namely, I'm not doing it by means of usual above declaration,)
      This is because I want to let a function create the pointer to the new objects,
      AggVertexes, and hand them back to this produce() method. Check this site:
      http://www.gotw.ca/publications/using_auto_ptr_effectively.htm. EC, 21-Nov-2010.
    */


    art::Handle< std::vector<recob::EndPoint2D> > epListHandle;
    evt.getByLabel(fEndPointModuleLabel,epListHandle);
    for(unsigned int ii = 0; ii < epListHandle->size(); ++ii)
      {
	art::Ptr<recob::EndPoint2D> ep(epListHandle, ii);
	feplist.push_back(ep); // class member
      }

    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrack3DModuleLabel,trackListHandle);
    for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii)
      {
	art::Ptr<recob::Track> track(trackListHandle, ii);
	ftracklist.push_back(track); // class member
      }


    // Only match strong vertices to tracks.
    art::PtrVector<recob::EndPoint2D>::const_iterator epIter = feplist.begin();

    while (epIter != feplist.end())  
      {            
	art::Ptr <recob::EndPoint2D> ep = (*epIter);  
	if (ep->ID() < 3) feplistStrong.push_back(ep); // -- cuz there's one less
	epIter++;
      }

    // We will suck up the hits out of each vertex and out of the tracks
    // and see if there's an overlap of a sufficient (>0) number of hits. If so,
    // call the track a match to the vtx. Then, .... stick the track pointer(s)
    // into the AggVertex object.  EC, 23-July-2010.

    std::auto_ptr< std::vector<recob::Vertex> > vcol (MatchV2T());

    evt.put(vcol);
  }

  std::auto_ptr< std::vector<recob::Vertex> >  AggregateVertex::MatchV2T()
  {
    mf::LogInfo("AggregateVertex") << "AggregateEvent::MatchV2T(): (strong) vertexlistStrong"
				   << " and tracklist lengths are " 
				   << feplistStrong.size() << " and " <<  ftracklist.size() << ".";


    // Bail if there are no tracks or vertices.
    //  if (!((int)vertexlistStrong.size()) || !((int)tracklist.size())) return NULL;
    if (feplistStrong.isNull() || ftracklist.isNull()) {
      return std::auto_ptr< std::vector<recob::Vertex> > (new std::vector<recob::Vertex>);
    }

    art::ServiceHandle<geo::Geometry> geom;
    // Loop on the vertices, and all the hits in each
    art::PtrVector<recob::EndPoint2D>::const_iterator epIter = feplistStrong.begin();
    std::auto_ptr< std::vector<recob::Vertex> > verts(new std::vector<recob::Vertex>);

    while (epIter != feplistStrong.end()){            
      art::PtrVector<recob::Track>  tlistAssoc; // Will fill with matching tracks.
      art::PtrVector<recob::Shower> slistAssoc; // Will fill with matching tracks.
      art::Ptr <recob::EndPoint2D> ep = (*epIter);
      art::PtrVector<recob::Hit> hitvertexlistStrong;
      
      art::PtrVector<recob::Hit> bv;
      for(unsigned int t = 0; t < geom->NTPC(); ++t){
	for(unsigned int ii = 0; ii < geom->Nplanes(t); ++ii){
	  geo::View_t view = geom->TPC(t).Plane(ii).View();
	  // std::cout << "AggregateEvent::MatchV2T(): Vtx hit view is " << view << std::endl;
	  bv = ep->Hits(view);
	  for (size_t jj=0; jj < bv.size(); jj++) hitvertexlistStrong.push_back(bv[jj]);
	}
      }// end loop over tpcs

      // Should be just one hit per vtx, as per Josh, but we loop anyway.
      art::PtrVector<recob::Hit>::const_iterator hvIter = hitvertexlistStrong.begin();
      while (hvIter != hitvertexlistStrong.end())  {
	art::Ptr<recob::Hit> hitv = (*hvIter);
	// Now loop on the tracks and all the clusters in each, and all
	// the hits in each of those.
	art::PtrVector<recob::Track>::const_iterator tIter = ftracklist.begin();
      
	while (tIter != ftracklist.end()){
	
	  art::Ptr<recob::Track> trk = (*tIter);
	
	  // Don't use this one, cuz these are Maddolena's
	  // added hits. I want the originals: the U,V hits.
	  // std::vector<const recob::Cluster*>clustlist3 = trk->Clusters(geo::k3D);
	  art::PtrVector<recob::Cluster>clustlistU = trk->Clusters(geo::kU);
	  art::PtrVector<recob::Cluster>clustlistV = trk->Clusters(geo::kV);
	  art::PtrVector<recob::Cluster>clustlist;
	  for (size_t jj=0; jj < clustlistU.size(); jj++) {clustlist.push_back(clustlistU[jj]);}
	  for (size_t jj=0; jj < clustlistV.size(); jj++) {clustlist.push_back(clustlistV[jj]);}
	
	  art::PtrVector<recob::Cluster>::const_iterator cIter = clustlist.begin();
	
	  while (cIter != clustlist.end()){
	    art::Ptr<recob::Cluster> cls = (*cIter);
	  
	    art::PtrVector<recob::Hit> hittlist;
	    art::PtrVector<recob::Hit> b;
	  
	    for(unsigned int t = 0; t < geom->NTPC(); ++t){
	      for(unsigned int ii = 0; ii < geom->TPC(t).Nplanes(); ++ii){
		geo::View_t view = geom->Plane(ii).View();
		b = cls->Hits(view);
		for (size_t jj=0; jj < b.size(); jj++) hittlist.push_back(b[jj]);
	      }
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
	  
	    art::PtrVector<recob::Hit>::const_iterator htIter = hittlist.begin();
	    while (htIter != hittlist.end())  {
	    
	      art::Ptr<recob::Hit> hitt = (*htIter);
	    
	      if (hitt==hitv){
		//std::cout << "AggregateEvent::MatchV2T(): BooYiggity! Vtx/Trk hit match! hitt, hitv= " 
		//<< hitt << " " <<hitv << std::endl;
		//std::cout << "AggregateEvent::MatchV2T(): vtx, trk " << vtx<< " " <<trk << std::endl;
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

      // Now if matching tracks were found for this vertex then create the recob::Vertex object.
    
      if (tlistAssoc.size()>0){
	// \todo Really need to also determine the xyz position of the found vertex
	// \todo Also should determine the ID for this vertex
	double xyz[3] = {-999., -999., -999.};
	verts->push_back(recob::Vertex(tlistAssoc, slistAssoc, xyz));
      }
	      
      epIter++;
      //HnhinVtxes->Fill(hitvertexlistStrong.size(),1);
	      
    }     // end epIter loop 
    

    return verts;
      
  } // end AggregateVertex::MatchV2T


}// end namespace
