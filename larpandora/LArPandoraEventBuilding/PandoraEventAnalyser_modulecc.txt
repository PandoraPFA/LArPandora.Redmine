////////////////////////////////////////////////////////////////////////
// Class:       PandoraEventAnalyser
// Plugin Type: analyzer (art v2_11_03)
// File:        PandoraEventAnalyser_module.cc
//
// Generated at Sun Oct 28 14:38:04 2018 by Wouter Van de pontseele using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

class PandoraEventAnalyser;


class PandoraEventAnalyser : public art::EDAnalyzer {
public:
  explicit PandoraEventAnalyser(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraEventAnalyser(PandoraEventAnalyser const &) = delete;
  PandoraEventAnalyser(PandoraEventAnalyser &&) = delete;
  PandoraEventAnalyser & operator = (PandoraEventAnalyser const &) = delete;
  PandoraEventAnalyser & operator = (PandoraEventAnalyser &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:
  int test;
  // Declare member data here.

};

DEFINE_ART_MODULE(PandoraEventAnalyser)

PandoraEventAnalyser::PandoraEventAnalyser(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void PandoraEventAnalyser::analyze(art::Event const & e)
{
  test = 1;
  std::cout << "Test: Hello World" << std::endl;
  // Implementation of required member function here.
}
