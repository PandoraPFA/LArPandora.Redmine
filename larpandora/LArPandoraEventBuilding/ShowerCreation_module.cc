/**
 *  @file   larpandora/LArPandoraEventBuilding/ShowerCreation_module.cc
 *
 *  @brief  module for shower creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace lar_pandora
{

class ShowerCreation : public art::EDProducer
{
public:
  explicit ShowerCreation(fhicl::ParameterSet const & p);

  ShowerCreation(ShowerCreation const &) = delete;
  ShowerCreation(ShowerCreation &&) = delete;
  ShowerCreation & operator = (ShowerCreation const &) = delete;
  ShowerCreation & operator = (ShowerCreation &&) = delete;

  void produce(art::Event & e) override;
};

DEFINE_ART_MODULE(ShowerCreation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "canvas/Utilities/InputTag.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace lar_pandora
{

ShowerCreation::ShowerCreation(fhicl::ParameterSet const & p)
{
    // Call appropriate produces<>() functions here.
}

void ShowerCreation::produce(art::Event & e)
{
    // Implementation of required member function here.
}

} // namespace lar_pandora
