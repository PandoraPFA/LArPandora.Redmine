////////////////////////////////////////////////////////////////////////
/// \file NUANCEGen.h
///
/// Module to produce neutrino interaction monte carlo using NUANCE EventFile
///
/// \version $Id: NUANCEGen.h,v 1.3 2010/02/15 19:10:40 brebel Exp $
/// \author  brebel@fnal.gov
/// \author  saima@ksu.edu
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_NUANCEGEN_H
#define EVGEN_NUANCEGEN_H

#include "art/Framework/Core/EDProducer.h"

#include "TStopwatch.h"

class TH1F;
class TH2F;

namespace simb { class MCTruth;     }

namespace evgen {
  /// A module to check the results from the Monte Carlo generator
  class NUANCEGen : public art::EDProducer {
  public:
    explicit NUANCEGen(fhicl::ParameterSet const &pset);
    virtual ~NUANCEGen();                        

    void produce(art::Event& evt);  
    void beginJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);
    void endJob();

  private:

        std::string ParticleStatus(int StatusCode);
        std::string ReactionChannel(int ccnc,int mode);
    
        void FillHistograms(simb::MCTruth mc);

	std::string         fNuanceFile;
	std::ifstream      *fEventFile;
	TStopwatch          fStopwatch;      ///keep track of how long it takes to run the job
	
	
	std::string fNUANCEModuleLabel;
	
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

#endif // EVGEN_NUANCEGEN_H
////////////////////////////////////////////////////////////////////////
