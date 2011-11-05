////////////////////////////////////////////////////////////////////////
/// \file CosmicsGen.h
/// 
/// Module to produce cosmic ray MC events using CRY
///
/// \version $Id: CosmicsGen.h,v 1.2 2010/02/15 19:10:40 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_COSMICSGEN_H
#define EVGEN_COSMICSGEN_H

#include "art/Framework/Core/EDProducer.h"

class TH1F;
class TH2F;

namespace evgb { class CRYHelper;   }  

namespace evgen {

  /// A module to check the results from the Monte Carlo generator
  class CosmicsGen : public art::EDProducer {
  public:
    explicit CosmicsGen(fhicl::ParameterSet const& pset);
    virtual ~CosmicsGen();                        


    void produce(art::Event& evt);  
    void beginJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    evgb::CRYHelper* fCRYHelp; ///< CRY generator object

    TH1F* fNhitHisto;          ///< Histogram of number of hits in spill	       
    TH1F* fDminHisto0;         ///< Closest approach for particles leaving no hits 
    TH1F* fDminHisto1;         ///< Closest approach for particles leaving hits    
    TH2F* fPhotonAngles;       ///< Photon rate vs angle			       
    TH2F* fPhotonAnglesLo;     ///< Photon rate vs angle, low momenta	       
    TH2F* fPhotonAnglesMi;     ///< Photon rate vs angle, middle momenta	       
    TH2F* fPhotonAnglesHi;     ///< Photon rate vs angle, high momenta	       
    TH1F* fPhotonCosQ;         ///< Photon rate vs cos(Q)			       
    TH1F* fPhotonEnergy;       ///< Photon energy (GeV)                            

    TH2F* fElectronAngles;     ///< Electron rate vs angle		 
    TH2F* fElectronAnglesLo;   ///< Electron rate vs angle, low momenta	 
    TH2F* fElectronAnglesMi;   ///< Electron rate vs angle, middle momenta 
    TH2F* fElectronAnglesHi;   ///< Electron rate vs angle, high momenta	 
    TH1F* fElectronCosQ;       ///< Electron rate vs cos(Q)		 
    TH1F* fElectronEnergy;     ///< Electron energy (GeV)                  

    TH2F* fMuonAngles;         ///< Muon rate vs angle		 
    TH2F* fMuonAnglesLo;       ///< Muon rate vs angle, low momenta	 
    TH2F* fMuonAnglesMi;       ///< Muon rate vs angle, middle momenta 
    TH2F* fMuonAnglesHi;       ///< Muon rate vs angle, high momenta	 
    TH1F* fMuonCosQ;           ///< Muon rate vs cos(Q)		 
    TH1F* fMuonEnergy;         ///< Muon energy (GeV)                  

    TH2F* fHitsYZ;             ///< hit locations in YZ view
    TH2F* fHitsXZ; 	       ///< hit locations in XZ view

  };
};

#endif // EVGEN_COSMICSGEN_H
////////////////////////////////////////////////////////////////////////
