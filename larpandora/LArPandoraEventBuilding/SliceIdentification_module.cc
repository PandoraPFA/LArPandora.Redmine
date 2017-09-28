/**
 *  @file   larpandora/LArPandoraEventBuilding/SliceIdentification_module.cc
 *
 *  @brief  module for slice identification and tagging
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace lar_pandora
{

class SliceIdentification : public art::EDProducer
{
public:
    explicit SliceIdentification(fhicl::ParameterSet const & p);

    SliceIdentification(SliceIdentification const &) = delete;
    SliceIdentification(SliceIdentification &&) = delete;
    SliceIdentification & operator = (SliceIdentification const &) = delete;
    SliceIdentification & operator = (SliceIdentification &&) = delete;

    void produce(art::Event & e) override;
};

DEFINE_ART_MODULE(SliceIdentification)

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

SliceIdentification::SliceIdentification(fhicl::ParameterSet const & p)
{
    // Call appropriate produces<>() functions here.
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceIdentification::produce(art::Event & e)
{
    // Implementation of required member function here.
}

} // namespace lar_pandora
