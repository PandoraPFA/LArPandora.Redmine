////////////////////////////////////////////////////////////////////////
/// \file GENIEGen.h
///
/// Module to produce neutrino interaction monte carlo using GENIE
///
/// \version $Id: GENIEGen.h,v 1.3 2010/02/15 19:10:40 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_GENIEGEN_H
#define EVGEN_GENIEGEN_H

#include "art/Framework/Core/EDProducer.h"

#include "TStopwatch.h"

class TH1F;
class TH2F;

namespace simb { class MCTruth;     }

namespace evgb { class GENIEHelper; }

///Event Generation using GENIE, cosmics or single particles
namespace evgen {
  /// A module to check the results from the Monte Carlo generator
  class GENIEGen : public art::EDProducer {
  public:
    explicit GENIEGen(fhicl::ParameterSet const &pset);
    virtual ~GENIEGen();                        

    void produce(art::Event& evt);  
    void beginJob();
    void beginRun(art::Run& run);
    void endSubRun(art::SubRun& sr);

  private:

    std::string ParticleStatus(int StatusCode);
    std::string ReactionChannel(int ccnc,int mode);
    
    void FillHistograms(simb::MCTruth mc);

    evgb::GENIEHelper  *fGENIEHelp;       ///< GENIEHelper object
    int                 fPassEmptySpills; ///< whether or not to kill evnets with no interactions
    TStopwatch          fStopwatch;       ///keep track of how long it takes to run the job

    TH1F* fGenerated[6];  ///< Spectra as generated

    TH1F* fVertexX;    ///< vertex location of generated events in x
    TH1F* fVertexY;    ///< vertex location of generated events in y
    TH1F* fVertexZ;    ///< vertex location of generated events in z

    TH2F* fVertexXY;   ///< vertex location in xy
    TH2F* fVertexXZ;   ///< vertex location in xz
    TH2F* fVertexYZ;   ///< vertex location in yz

    TH1F* fDCosX;      ///< direction cosine in x
    TH1F* fDCosY;      ///< direction cosine in y
    TH1F* fDCosZ;      ///< direction cosine in z

    TH1F* fMuMomentum; ///< momentum of outgoing muons
    TH1F* fMuDCosX;    ///< direction cosine of outgoing mu in x
    TH1F* fMuDCosY;    ///< direction cosine of outgoing mu in y
    TH1F* fMuDCosZ;    ///< direction cosine of outgoing mu in z

    TH1F* fEMomentum;  ///< momentum of outgoing electrons
    TH1F* fEDCosX;     ///< direction cosine of outgoing e in x
    TH1F* fEDCosY;     ///< direction cosine of outgoing e in y
    TH1F* fEDCosZ;     ///< direction cosine of outgoing e in z

    TH1F* fCCMode;      ///< CC interaction mode
    TH1F* fNCMode;      ///< CC interaction mode

    TH1F* fDeltaE;     ///< difference in neutrino energy from MCTruth::Enu() vs TParticle
    TH1F* fECons;      ///< histogram to determine if energy is conserved in the event

  };
};

#endif // MCCHK_GENCHECK_H
////////////////////////////////////////////////////////////////////////
