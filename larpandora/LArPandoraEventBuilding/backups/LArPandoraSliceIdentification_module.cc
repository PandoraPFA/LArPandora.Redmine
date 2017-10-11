/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraSliceIdentification_module.cc
 *
 *  @brief  module for lar pandora slice identification and tagging
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraSliceIdentification : public art::EDProducer
{
public:
    explicit LArPandoraSliceIdentification(fhicl::ParameterSet const & p);

    LArPandoraSliceIdentification(LArPandoraSliceIdentification const &) = delete;
    LArPandoraSliceIdentification(LArPandoraSliceIdentification &&) = delete;
    LArPandoraSliceIdentification & operator = (LArPandoraSliceIdentification const &) = delete;
    LArPandoraSliceIdentification & operator = (LArPandoraSliceIdentification &&) = delete;

    void produce(art::Event & e) override;
};

DEFINE_ART_MODULE(LArPandoraSliceIdentification)

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

LArPandoraSliceIdentification::LArPandoraSliceIdentification(fhicl::ParameterSet const & p)
{
    // Call appropriate produces<>() functions here.
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraSliceIdentification::produce(art::Event & e)
{
    // Implementation of required member function here.
}

} // namespace lar_pandora
