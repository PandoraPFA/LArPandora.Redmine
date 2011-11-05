#ifndef EventMaker_h
#define EventMaker_h
////////////////////////////////////////////////////////////////////////
// Class:       EventMaker
// Module Type: producer
// File:        EventMaker.h
//
// Generated at Wed Jul 13 15:06:57 2011 by Brian Rebel using artmod
// from art v0_07_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"

namespace event {
  class EventMaker;
}

class event::EventMaker : public art::EDProducer {
public:
  explicit EventMaker(fhicl::ParameterSet const &p);
  virtual ~EventMaker();

  virtual void produce(art::Event &e);

  virtual void reconfigure(fhicl::ParameterSet const & p);

private:

  std::string fVertexModuleLabel;  ///< label of the module making the recob::Vertex objects
  double      fProximity;          ///< how close a vertex needs to be to another to be from the same event

};
#endif /* EventMaker_h */
