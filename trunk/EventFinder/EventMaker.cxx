////////////////////////////////////////////////////////////////////////
// Class:       EventMaker
// Module Type: producer
// File:        EventMaker.h
//
// Generated at Wed Jul 13 15:06:57 2011 by Brian Rebel using artmod
// from art v0_07_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Principal/Event.h"
#include "EventFinder/EventMaker.h"
#include "RecoBase/recobase.h"

//--------------------------------------------------------
event::EventMaker::EventMaker(fhicl::ParameterSet const &p)
{
  this->reconfigure(p);
  
  produces< std::vector<recob::Event> >();
}

//--------------------------------------------------------
event::EventMaker::~EventMaker() 
{
  // Clean up dynamic memory and other resources here.
}

//--------------------------------------------------------
void event::EventMaker::produce(art::Event &e) 
{

  // make the auto_ptr of the vector for the recob::Events
  std::auto_ptr< std::vector<recob::Event> > eventcol(new std::vector<recob::Event>);

  // first get the recob::Vertex objects out of the event
  art::Handle< std::vector<recob::Vertex> > vtxHandle;
  e.getByLabel(fVertexModuleLabel, vtxHandle);

  // get a collection of art::Ptrs
  std::list< art::Ptr<recob::Vertex> > vtxs;
  art::fill_ptr_list(vtxs, vtxHandle);

  // if only 1 vertex in the event, life is easy
  if(vtxs.size() == 1){
    art::PtrVector<recob::Vertex> vtxvec;
    vtxvec.push_back(*(vtxs.begin()));
    recob::Event evt(vtxvec, 0);
    eventcol->push_back(evt);
    e.put(eventcol);
    return;
  }

  // now the hard part, we have multiple vertex objects
  // things to consider:
  // 1. all particles coming from a common vertex location are in the same event
  // 2. particles coming from decays of another particle, ie the electron from a muon decay
  // 3. pi-zero decay vertex will be offset from the interaction vertex
  int evtctr = 0;
  std::list< art::Ptr<recob::Vertex> >::iterator itr  = vtxs.begin();
  std::list< art::Ptr<recob::Vertex> >::iterator itr2 = vtxs.begin();
  while( itr != vtxs.end() ){
    
    art::PtrVector< recob::Vertex > vtxVec;

    // get the current vertex object and put it into the vector
    art::Ptr<recob::Vertex> curvtx = *itr;
    vtxVec.push_back(curvtx);

    double curxyz[3] = {0.};
    curvtx->XYZ(curxyz);

    // make itr2 the next one in the list, also remove this one from the list
    // as it is being used
    itr2 = vtxs.erase(itr);
    while( itr2 != vtxs.end() ){
      art::Ptr<recob::Vertex> nexvtx = *itr2;
      
      // get the xyz location of the vertex to compare to the current vertex
      double nextxyz[3] = {0.};
      nexvtx->XYZ(nextxyz);
      if( sqrt((curxyz[0]-nextxyz[0])*(curxyz[0]-nextxyz[0]) +
	       (curxyz[1]-nextxyz[1])*(curxyz[1]-nextxyz[1]) +
	       (curxyz[2]-nextxyz[2])*(curxyz[2]-nextxyz[2]) ) <= fProximity){
	
	// add this one to the vector and remove it from the list as it is being used
	vtxVec.push_back(nexvtx);
	itr2 = vtxs.erase(itr2);
      }// end if vertices are close enough to be considered from the same event
      else itr2++;
    }// end inner loop
    
    // make an event from these vertex objects and add them to the collection
    recob::Event evt(vtxVec, ++evtctr);
    eventcol->push_back(evt);

    // move the initial iterator forward
    itr++;

  }

  // put the collection of events in the art::Event
  e.put(eventcol);

  return;

}

//--------------------------------------------------------
void event::EventMaker::reconfigure(fhicl::ParameterSet const & p) 
{
  fVertexModuleLabel = p.get< std::string >("VertexModuleLabel");
  fProximity         = p.get< double      >("Proximity");
}

