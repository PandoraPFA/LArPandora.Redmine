///////////////////////////////////////////////////////////////////////
/// \file    VertexCheater.h
/// \brief   make vertices using MC truth information
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef VERTEX_VERTEXCHEATER_H
#define VERTEX_VERTEXCHEATER_H
#include <string>

#include "art/Framework/Core/EDProducer.h"

namespace vertex {
  class VertexCheater : public art::EDProducer {
  public:
    explicit VertexCheater(fhicl::ParameterSet const& pset);
    virtual ~VertexCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fCheatedTrackLabel;  ///< label for module creating recob::Prong objects
    std::string fCheatedShowerLabel; ///< label for module creating recob::Prong objects
    std::string fG4ModuleLabel;      ///< label for module running G4 and making particles, etc

  };
}
#endif
