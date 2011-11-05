////////////////////////////////////////////////////////////////////////
/// \file FileMuons.h
///
/// Module to produce a set list of particles for a MC event
///
/// \version $Id: FileMuons.h,v 1.2 2010/02/15 19:10:40 echurch Exp $
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
// Framework includes.
#include <vector>
#include <string>
#include "art/Framework/Core/EDProducer.h"
#include "TFile.h"
#include "TTree.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

namespace simb { class MCTruth; }

namespace evgen {

  /// module to produce single or multiple specified particles in the detector
  class FileMuons : public art::EDProducer {

  public:
    explicit FileMuons(fhicl::ParameterSet const& pset);
    virtual ~FileMuons();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginJob();
    void beginRun(art::Run& run);
    void endJob();

  private:

    void ReadEvents(simb::MCTruth &mct);        

    int                 fEventNumberOffset;  // Where in file to start.
    int                 fSeed;           // random number seed    
    std::vector<int>    fPDG;           
    std::string fFileName;
    std::string fMuonsFileType;
    std::string fTreeName;
    std::vector<std::string> fBranchNames;

    ifstream *fMuonFile;
    TFile *fMuonFileR;
    TTree *TNtuple;
    unsigned int countFile;

    Float_t xtmp, ytmp, ztmp;
    Float_t pxtmp, pytmp, pztmp;
    Float_t charge;
    Float_t         E;
    Float_t         costheta;
    Float_t         phi;
    Float_t         xdet;
    Float_t         ydet;
    Float_t         zdet;

    TBranch        *b_x;   //!
    TBranch        *b_y;   //!
    TBranch        *b_z;   //!
    TBranch        *b_E;   //!
    TBranch        *b_costheta;   //!
    TBranch        *b_phi;   //!
    TBranch        *b_xdet;   //!
    TBranch        *b_ydet;   //!
    TBranch        *b_zdet;   //!
    TBranch        *b_px;   //!
    TBranch        *b_py;   //!
    TBranch        *b_pz;   //!
    TBranch        *b_charge;   //!


  };
};
////////////////////////////////////////////////////////////////////////
