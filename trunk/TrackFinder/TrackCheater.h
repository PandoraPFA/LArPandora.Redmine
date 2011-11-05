///////////////////////////////////////////////////////////////////////
/// \file    TrackCheater.h
/// \brief   make prongs using MC truth information
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef TRKF_TRACKCHEATER_H
#define TRKF_TRACKCHEATER_H
#include <string>

#include "art/Framework/Core/EDProducer.h"

namespace trkf {
  class TrackCheater : public art::EDProducer {
  public:
    explicit TrackCheater(fhicl::ParameterSet const& pset);
    virtual ~TrackCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fCheatedClusterLabel; ///< label for module creating recob::Cluster objects	   
    std::string fG4ModuleLabel;       ///< label for module running G4 and making particles, etc

  };
}
#endif
