///////////////////////////////////////////////////////////////////////
/// \file    ShowerCheater.h
/// \brief   make prongs using MC truth information
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef SHWF_SHOWERCHEATER_H
#define SHWF_SHOWERCHEATER_H
#include <string>

#include "art/Framework/Core/EDProducer.h"

namespace shwf {
  class ShowerCheater : public art::EDProducer {
  public:
    explicit ShowerCheater(fhicl::ParameterSet const& pset);
    virtual ~ShowerCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fCheatedClusterLabel; ///< label for module creating recob::Cluster objects	   
    std::string fG4ModuleLabel;       ///< label for module running G4 and making particles, etc

  };
}
#endif
