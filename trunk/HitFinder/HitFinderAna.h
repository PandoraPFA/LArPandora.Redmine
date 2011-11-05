////////////////////////////////////////////////////////////////////////
// $Id: HitFinder.h,v 1.12 2010/04/23 20:30:53 seligman Exp $
//
// HitFinder class designed to simulate signal on a wire in the TPC
//
// echurch@fnal.gov, brebel@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////
#ifndef HITFINDERANA_H
#define HITFINDERANA_H


#include "RecoBase/recobase.h"
#include "Utilities/LArFFT.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires. 
  class HitFinderAna : public art::EDAnalyzer {
    
  public:
        
    explicit HitFinderAna(fhicl::ParameterSet const& pset); 
    virtual ~HitFinderAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::string            fFFTHitFinderModuleLabel;
    std::string            fLArG4ModuleLabel;
    
      TTree* fHTree;
      Int_t fRun;
      Int_t fEvt;
      Int_t fNp0;
      Int_t fNp1;
      Int_t fNp2;
      Int_t fN3p0;
      Int_t fN3p1;
      Int_t fN3p2;
      Float_t* fTimep0;
      Float_t* fTimep1;
      Float_t* fTimep2;
      Int_t* fWirep0;
      Int_t* fWirep1;
      Int_t* fWirep2;
      Float_t* fChgp0;
      Float_t* fChgp1;
      Float_t* fChgp2;
      Float_t* fXYZp0;
      Float_t* fXYZp1;
      Float_t* fXYZp2;

      Int_t*  fMCPdg0;
      Int_t*  fMCTId0;
      Float_t*  fMCE0;
      Int_t*  fMCPdg1;
      Int_t*  fMCTId1;
      Float_t*  fMCE1;
      Int_t*  fMCPdg2;
      Int_t*  fMCTId2;
      Float_t*  fMCE2;

  }; // class HitFinderAna

}

#endif // HITFINDERANA_H
