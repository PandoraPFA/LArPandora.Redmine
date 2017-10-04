////////////////////////////////////////////////////////////////////////
// Class:       LArPandoraNeutrinoId
// Plugin Type: producer (art v2_07_03)
// File:        LArPandoraNeutrinoId_module.cc
//
// Generated at Mon Oct  2 11:54:57 2017 by Andrew D. Smith using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraSlices.h"

#include <memory>

class LArPandoraNeutrinoId;


class LArPandoraNeutrinoId : public art::EDProducer {
public:
  explicit LArPandoraNeutrinoId(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LArPandoraNeutrinoId(LArPandoraNeutrinoId const &) = delete;
  LArPandoraNeutrinoId(LArPandoraNeutrinoId &&) = delete;
  LArPandoraNeutrinoId & operator = (LArPandoraNeutrinoId const &) = delete;
  LArPandoraNeutrinoId & operator = (LArPandoraNeutrinoId &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.

};


LArPandoraNeutrinoId::LArPandoraNeutrinoId(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

void LArPandoraNeutrinoId::produce(art::Event & e)
{
  // Implementation of required member function here.
}

void LArPandoraNeutrinoId::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(LArPandoraNeutrinoId)
