////////////////////////////////////////////////////////////////////////
/// file LambdaAna.h
//
/// author saima@ksu.edu
/// 
////////////////////////////////////////////////////////////////////////
#ifndef CCHYP_LAMBDAANA_H
#define CCHYP_LAMBDAANA_H
//#include "JobControl/Module.h"
#include "art/Framework/Core/EDAnalyzer.h"

//namespace edm { class EventHandle; }
//namespace cfg { class Config;      }
class TH1F;
class TH2F;
class TGraph;
class TFile;

namespace cchyp {
  /// A module to check the results from the Monte Carlo generator
  class LambdaAna : public art::EDAnalyzer {
  public:
    explicit LambdaAna(fhicl::ParameterSet const& pset);
    virtual ~LambdaAna();                        
    //void Update(const cfg::Config& c);
    void analyze(const art::Event&);
    void beginJob(); 
    // void EndJob(const edm::Event&, const edm::EventSetup &); 
    // jobc::Result Ana(const edm::EventHandle& evt);
  
  private:
    
    std::string fGenieModuleLabel;
    std::string fLArG4ModuleLabel; 
    std::string fHitModuleLabel;
    std::string fDetSimModuleLabel;
    std::string fLineMergerModuleLabel; 
    std::string fTrackModuleLabel; 
    std::string fVertexModuleLabel;

    TH1F* fXreco_Xtrue;
    TH1F* fYreco_Ytrue;
    TH1F* fZreco_Ztrue;

    TH1F* fXreco_Xmuon;
    TH1F* fYreco_Ymuon;
    TH1F* fZreco_Zmuon;

    TH1F* fMINDIST;
    TH2F* fRtracksEnu;
    TH2F* fTtracksEnu;
    TH2F* fRtracksEp;
    TH2F* fTtracksEp;

    TH2F* fTRtracks;

    TH1F* fRVertexX;   
    TH1F* fRVertexY;   
    TH1F* fRVertexZ;   

    TH2F* fRVertexXY;  
    TH2F* fRVertexXZ;  
    TH2F* fRVertexYZ; 

    TH1F* fRmuonX_TmuonX; 
    TH1F* fRmuonY_TmuonY; 
    TH1F* fRmuonZ_TmuonZ; 
    
    TH1F*   fTVertexX;
    TH1F*   fTVertexY;
    TH1F*   fTVertexZ; 

    TH2F* fXtrue_vs_Xreco;
    TH2F* fYtrue_vs_Yreco;
    TH2F* fZtrue_vs_Zreco;

    TH1F* fEventsWith_greaterTHAN_2cmX; 
    TH1F* fEventsWith_greaterTHAN_2cmY; 
    TH1F* fEventsWith_greaterTHAN_2cmZ; 

    TH2F* fXtrue_vs_Xmuon;
    TH2F* fYtrue_vs_Ymuon;
    TH2F* fZtrue_vs_Zmuon;
    
    TH2F* fXmuon_vs_Xreco;
    TH2F* fYmuon_vs_Yreco;
    TH2F* fZmuon_vs_Zreco;
    TH1F* fnumber_true; 

    TH1F* fTotal_CCQE;
    TH1F* fCCQE_FV; 
    TH1F* fCCQE_decay_FV; 
    TH1F* fMuon_Escapes; 
    TH1F* fMuon_Enters_MINOS; 
    TH1F* fProton_Escapes; 
    TH1F* fPion_Escapes;

    TH1F* fMuon_Length; 
    TH1F* fProton_Length;
    TH1F* fPion_Length; 
      
    TH2F* fMuon_KE_Length; 
    TH2F* fProton_KE_Length;
    TH2F* fPion_KE_Length;

    TH1F* fU_no_clus;
    TH1F* fV_no_clus;

    TH1F* fU_clu_length; 
    TH1F* fV_clu_length;
    TH1F* fUwire_span_length;
    TH1F* fVwire_span_length;

    TH2F* fMuon_Reco_Eff; 
    TH2F* fProtonn_Reco_Eff;
    TH2F* fPion_Reco_Eff;

    TH2F* fMuon_Reco_Eff_vs_hits;
    TH2F* fProtonn_Reco_Eff_vs_hits;
    TH2F* fPion_Reco_Eff_vs_hits;

    TH1F* flambda_decay_length;
    TH1F* fmuon_lambda_angle;
    TH1F* fproton_pion_angle;
    TH1F* flambda_KE;

    TH2F* fMuon_Reco_Eff_3d; 
    TH2F* fProtonn_Reco_Eff_3d;
    TH2F* fPion_Reco_Eff_3d;
    
    TH1F* fNumber_tracks;

    TH2F* fpurity_muon;
    TH2F* fpurity_proton;
    TH2F* fpurity_pion;

    TH2F* fcompleteness_muon;
    TH2F* fcompleteness_proton;
    TH2F* fcompleteness_pion;

     TH2F* fmuon_cluster_hits_UV;
     TH2F* fproton_cluster_hits_UV;
     TH2F* fpion_cluster_hits_UV;

     TH1F* fmuon_reco_percentage;
     TH1F* fproton_reco_percentage;
     TH1F* fpion_reco_percentage;

     TH2F* fmuon_cluster_charge_UV;
     TH2F* fproton_cluster_charge_UV;
     TH2F* fpion_cluster_charge_UV;
      
  };
};

#endif // CCHYP_LAMBDAANA_H
////////////////////////////////////////////////////////////////////////
