////////////////////////////////////////////////////////////////////////
//
//   \file Track3DKalman.h
//
//   soderber@fnal.gov
//   kinga.partyka@yale.edu
//   joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef TRACK3DRECO_H
#define TRACK3DRECO_H

#include "art/Framework/Core/EDProducer.h"
#include <TTree.h>
#include <TMatrixT.h>

#include "Genfit/GFAbsTrackRep.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <vector>
#include <string>

//#include "RecoBase/SpacePoint.h"

namespace trkf {

  class Track3DKalman : public art::EDProducer {
    
  public:
    
    explicit Track3DKalman(fhicl::ParameterSet const& pset);
    ~Track3DKalman();
    
    //////////////////////////////////////////////////////////
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:
        
    std::string     fSpacePtsModuleLabel;// label for input collection
    std::string     fGenieGenModuleLabel;// label for input MC single particle generator
    std::string     fG4ModuleLabel;// label for input MC single particle generator
    bool fGenfPRINT;
      
    TFile *fileGENFIT;
    TTree *tree;

    TMatrixT<Double_t> *stMCT;
    TMatrixT<Double_t> *covMCT;
    TMatrixT<Double_t> *stREC;
    TMatrixT<Double_t> *covREC;
    Float_t chi2;
    Float_t chi2ndf;
    
    Float_t *fpRECt3D;    
    Float_t *fpRECL;
    Float_t *fpREC;
    Float_t *fpMCT;
    int nfail;
    int ndf;
    unsigned int evtt;
    unsigned int nTrks;
    unsigned int fptsNo;
    Float_t *fshx;
    Float_t *fshy;
    Float_t *fshz;
    unsigned int fDimSize; // if necessary will get this from pset in constructor.
  
    std::vector<double> fPosErr;
    std::vector<double> fMomErr;
    std::vector<double> fMomStart;
    genf::GFAbsTrackRep *repMC;
    genf::GFAbsTrackRep *rep;

  protected: 
    
  
  }; // class Track3DKalman

}

#endif // TRACK3DRECO_H
