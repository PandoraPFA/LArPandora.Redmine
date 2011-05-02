////////////////////////////////////////////////////////////////////////
//
// 
// \author tjyang@fnal.gov
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TTree;

///Track finding and building 
namespace t962 {

   
  class AnalysisTree : public art::EDAnalyzer {

  public:
          
    explicit AnalysisTree(fhicl::ParameterSet const& pset); 
    virtual ~AnalysisTree();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

  private:
    
    void ResetVars();

    TTree* fTree;
    //run information
    int run;
    int event;
    int isdata;
    //reconstructed information
    double vtxx;
    double vtxy;
    double vtxz;
    int ntrks;
    int nclusu;
    int nclusv;
    int nclusw;
    //matching information
    int matched;        //number of matched tracks
    double mtrk_mom;
    double mtrk_charge;
    double mtrk_dcosx;
    double mtrk_dcosy;
    double mtrk_dcosz;
    //mctruth information
    int inu;
    int ccnc;
    int mode;
    double enu;
    double Q2;
    double W;
    double nuvtxx;
    double nuvtxy;
    double nuvtxz;
    double lep_mom;
    double lep_dcosx;
    double lep_dcosy;
    double lep_dcosz;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fClusterModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fVertexModuleLabel;
    std::string fMINOSModuleLabel;
    std::string fTrackMatchModuleLabel;

  }; // class AnalysisTree

}

#endif 
