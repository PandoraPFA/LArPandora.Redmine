///////////////////////////////////////////////////////////////////////
/// \file    ClusterCheater.h
/// \brief   make clusters using MC truth information
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef CLUSTER_CLUSTERCHEATER_H
#define CLUSTER_CLUSTERCHEATER_H
#include <string>

#include "art/Framework/Core/EDProducer.h"

namespace cluster {
  class ClusterCheater : public art::EDProducer {
  public:
    explicit ClusterCheater(fhicl::ParameterSet const& pset);
    virtual ~ClusterCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fHitModuleLabel;    ///< label for module creating recob::Hit objects	   
    std::string fG4ModuleLabel;     ///< label for module running G4 and making particles, etc

  };
}
#endif
