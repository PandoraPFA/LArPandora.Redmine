////////////////////////////////////////////////////////////////////////
/// \file  ChargedPion.h
/// \module to extract PDG codes, energy, etc from mc events with pi+
///
/// \version $Id: LArG4.h,v 1.11 2010/06/04 21:47:27 bjpjones Exp $?
/// \author  ellen.klein@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef PiPlus_PiPlus_H
#define PiPlus_PiPlus_H 

#include "art/Framework/Core/EDAnalyzer.h"
#include <cstring>
#include <TTree.h>

namespace cpion {

  class ChargedPion : public art::EDAnalyzer{
  public:
    /// Standard constructor and destructor for an FMWK module.
    explicit ChargedPion(fhicl::ParameterSet const& pset);
    virtual ~ChargedPion();

    void analyze (const art::Event& evt); 
    void beginJob();

  private:

    std::string fG4ModuleLabel;  // name of the process module label that produced the input simb::MCTruths
    std::string fGenieModuleLabel;

    TTree *nTree;
    TTree *gTree;
    TTree *tTree;

    Int_t Evt; //Event number, used for all trees
    
    //Stuff i need for trees...
    Int_t pTNds;
    Int_t pTNdsOriginal;

    //Neutrino information, used for nTree. Pulled from Genie
    Int_t nPdg;
    Int_t nType;
    Double_t nVx; 
    Double_t nVy; 
    Double_t nVz; 
    Double_t nEnergy;   
    Double_t nQSquared;

    //Genie particle information, used for gTree
    Int_t gPdg;
    Int_t gTrackId;
    Int_t gMother;
    Int_t gStatus;
    Double_t gEnergy;
    Double_t gVx;
    Double_t gVy;
    Double_t gVz;


    //Geant Pion information, used for tTree
    Int_t tPdg;
    Int_t tTrackId;
    Int_t tMother;
    Int_t tDaughters;
    Int_t tDecay;
    Int_t tInelastic;
    Int_t tCapture;
    Double_t tEnergy;
    Double_t tVx;
    Double_t tVy;
    Double_t tVz;
    Double_t tEndx;
    Double_t tEndy;
    Double_t tEndz;
    Double_t tPx;
    Double_t tPy;
    Double_t tPz;
    Double_t tP;
    Double_t tPt;
    Double_t tDaughPdg[200];
    Double_t tProtonE[10];
    Double_t tNeutronE[10];
    Double_t tProtonDist[10];
    Int_t tProtons;
    Int_t tNeutrons;
    Double_t tTotalEnergy;
    Double_t tChargedEnergy;
    Double_t tTotalPPar;
    Double_t tTotalPPerp;
    Double_t tChargedPPar;
    Double_t tChargedPPerp;
    Int_t tNumPions;
    //-------
    Int_t tMuPdg;
    Double_t tMuEnergy;
    Double_t tMuVx;
    Double_t tMuVy;
    Double_t tMuVz;
    Double_t tMuEndx;
    Double_t tMuEndy;
    Double_t tMuEndz;
    Double_t tMuPx;
    Double_t tMuPy;
    Double_t tMuPz;
    Int_t tnumPrimary;
    Int_t numPiPlus; //NOT DECAY// over all run, not in tree just cout'd
    Int_t numPiMinus; //NOT DECAY// over all run, not in tree just cout'd
    Int_t num2Nu;
    Int_t LastEvt;
    Int_t PiLock;
  };

} // namespace PiPlus

#endif // PiPlus_PiPlus_H
//std::cout<<"TETING"<<std::endl;
