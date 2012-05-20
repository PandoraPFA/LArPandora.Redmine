////////////////////////////////////////////////////////////////////////
/// \file  TGMuon.h
/// \brief Generator for through-going muons
/// Module designed to produce through-going muons based on MINOS data for ArgoNeuT (run 3, neutrino-mode right now)
///
/// \author  joshua.spitz@yale.edu

#include <vector>
#include <string>
#include "art/Framework/Core/EDProducer.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"

namespace simb { class MCTruth; }

namespace t962 {

  /// module to produce single or multiple specified particles in the detector
  class TGMuon : public art::EDProducer {

  public:
    explicit TGMuon(fhicl::ParameterSet const& pset);
    virtual ~TGMuon();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void endSubRun(art::SubRun& sr);
    void reconfigure(fhicl::ParameterSet const& pset) {};

  private:

    void Sample(simb::MCTruth &mct);        
    
    TRandom3           *fRand;           ///< random number
    int                 fSeed;           ///< random number seed    
    int frun;
    int fsubRun;
    int fsnarl;
    int fntrkstp;
    float ftrkIndex;
    float ftrtgtd;
    float ftrkVtxX;
    float ftrkVtxY;
    float ftrkVtxZ;
    float ftrkdcosx;
    float ftrkdcosy;
    float ftrkdcosz; 
    float ftrkmom;
    float ftrkqp;
    float totalpot;

  };
};
////////////////////////////////////////////////////////////////////////
