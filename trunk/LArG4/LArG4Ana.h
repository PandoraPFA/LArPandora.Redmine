////////////////////////////////////////////////////////////////////////
/// \file  LArG4Ana.h
/// \brief Check of Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.h,v 1.11 2010/06/04 21:47:27 bjpjones Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_LARG4ANA_H
#define LARG4_LARG4ANA_H 

#include "art/Framework/Core/EDAnalyzer.h"

#include <cstring>
#include <TTree.h>
class TH1D;
class TH2D;

namespace simb{
  class MCTruth;
}

namespace sim{
  class LArVoxelList;
  class ParticleList;
}

///Geant4 interface 
namespace larg4 {  
 
  class LArG4Ana : public art::EDAnalyzer{
  public:
 
    /// Standard constructor and destructor for an FMWK module.
    explicit LArG4Ana(fhicl::ParameterSet const& pset);
    virtual ~LArG4Ana();

    void analyze (const art::Event& evt); 
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& pset);

  private:

    std::string fG4ModuleLabel;     ///< module label for the Geant
    std::string fTruthModuleLabel;  ///< module label for the Geant

    void CheckContainment(sim::LArVoxelList const& vlist,
			  const simb::MCTruth*     mct);

    TH1D *fDistFromVertex;     ///< distance of energy depositions from interaction vertex
    TH2D *fEnergyWithin5cm;    ///< fraction of energy deposited within 5  cm of shower axis vs total energy
    TH2D *fEnergyWithin10cm;   ///< fraction of energy deposited within 10 cm of shower axis vs total energy
    TH2D *fEnergyWithin15cm;   ///< fraction of energy deposited within 15 cm of shower axis vs total energy
    TH2D *fEnergyWithin20cm;   ///< fraction of energy deposited within 20 cm of shower axis vs total energy
    TH2D *fEnergyWithin25cm;   ///< fraction of energy deposited within 25 cm of shower axis vs total energy
    TH2D *fEnergyWithin30cm;   ///< fraction of energy deposited within 30 cm of shower axis vs total energy
    TH2D *fEnergyWithin35cm;   ///< fraction of energy deposited within 35 cm of shower axis vs total energy
    TH2D *fEnergyWithin40cm;   ///< fraction of energy deposited within 40 cm of shower axis vs total energy
    TH2D *fEnergyWithin45cm;   ///< fraction of energy deposited within 45 cm of shower axis vs total energy
    TH2D *fEnergyWithin50cm;   ///< fraction of energy deposited within 50 cm of shower axis vs total energy

    TH1D *fPDGCodes;
    TH1D *fPi0Momentum;
    TH1D *fnEnergy;
    TH1D *fnDist;
    
    //    Int_t stringDim = 35;

    TTree *fTree;
    Int_t fTEvt;
    Int_t fTSub;
    Int_t fTRun;
    Int_t fTPdg;
    Int_t fTID;
    Int_t fTNds;
    Int_t fTNdsOriginal;
    Int_t fTNds4;
    Int_t *fTDID;
    Int_t *fTDPdg;
    Float_t *fTDWt;
    Char_t fTProcess[35];
    Char_t fTVolume[35];
    Char_t fTTVolume[35]; // Termination Volume
    Char_t fTMaterial[35];
    Char_t fTDProcess[200][35];
    Int_t fTParentID; 
    Int_t fTStatus;
    Float_t fTWeight;
    Float_t* fT4Origin;
    Float_t* fT4DOrigin;
    Float_t* fT4Termination; // Termination Coordinates
    Float_t* fT4Momentum;
    Float_t* fT4DMomentum;
  };

} // namespace larg4

#endif // LARG4_LARG4_H
