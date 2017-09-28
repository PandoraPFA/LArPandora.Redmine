/**
 *  @file   larpandora/LArPandoraEventBuilding/TrackCreation_module.cc
 *
 *  @brief  module for track creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace lar_pandora
{

class TrackCreation : public art::EDProducer
{
public:
  explicit TrackCreation(fhicl::ParameterSet const & p);

  TrackCreation(TrackCreation const &) = delete;
  TrackCreation(TrackCreation &&) = delete;
  TrackCreation & operator = (TrackCreation const &) = delete;
  TrackCreation & operator = (TrackCreation &&) = delete;

  void produce(art::Event & e) override;
};

DEFINE_ART_MODULE(TrackCreation)

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

TrackCreation::TrackCreation(fhicl::ParameterSet const & p)
{
    // Call appropriate produces<>() functions here.
}

void TrackCreation::produce(art::Event & e)
{
    // Implementation of required member function here.
}

} // namespace lar_pandora
