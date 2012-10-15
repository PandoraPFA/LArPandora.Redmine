////////////////////////////////////////////////////////////////////////
//
// MuonAna class
//
// msoderbe@syr.edu
//
//  
////////////////////////////////////////////////////////////////////////

#ifndef MUONANA_H
#define MUONANA_H
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Core/FindMany.h"

#include "RecoBase/Track.h"
#include "RecoBase/Hit.h"
#include "T962/T962_Objects/MINOS.h"

#include <iostream>
#include <vector>
#include <string>

#include <TH1F.h>
#include <TH2F.h>
#include "TMath.h"


///T962 muon analysis code
namespace muons {
   
  class MuonAna :  public art::EDAnalyzer {
    
  public:
    
    explicit MuonAna(fhicl::ParameterSet const& pset); 
    virtual ~MuonAna();        

    void analyze (const art::Event& evt);
    void beginJob();

    bool BeginsOnBoundary(art::Ptr<recob::Track> lar_track);
    bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);
    void Compare(art::Ptr<recob::Track> lar_track, t962::MINOS minos_track,   
                 double &dx, double &dy, double &rdiff);
        
  private:
       
    //plots for TG Muons
    TH1F*  fDirX_start;
    TH1F*  fDirY_start;
    TH1F*  fDirZ_start;
    TH1F*  fCosX_start;
    TH1F*  fCosY_start;
    TH1F*  fCosZ_start;
    TH1F*  fX_start;
    TH1F*  fY_start;
    TH1F*  fZ_start;
    TH1F*  fX_end;
    TH1F*  fY_end;
    TH1F*  fZ_end;
    TH1F*  fTrackLength;
    TH1F*  fTheta;
    TH1F*  fPhi;
    TH1F*  fVertAngle;
    TH1F*  fHorizAngle;
    TH2F*  fStartXvsStartY;
    TH2F*  fStartZvsStartX;
    TH2F*  fStartZvsStartY;
    TH2F*  fEndXvsEndY;
    TH2F*  fEndZvsEndX;
    TH2F*  fEndZvsEndY;

    TH2F*  fChannelVsHitAmplitude;
    //TH2F*  fChannelVsHitAmplitude_Corrected;

    //plots for TG Muons matched to a pos/neg MINOS track
    TH1F*  fDirX_start_Pos;
    TH1F*  fDirY_start_Pos;
    TH1F*  fDirZ_start_Pos;
    TH1F*  fCosX_start_Pos;
    TH1F*  fCosY_start_Pos;
    TH1F*  fCosZ_start_Pos;
    TH1F*  fX_start_Pos;
    TH1F*  fY_start_Pos;
    TH1F*  fZ_start_Pos;
    TH1F*  fX_end_Pos;
    TH1F*  fY_end_Pos;
    TH1F*  fZ_end_Pos;
    TH1F*  fTrackLength_Pos;
    TH1F*  fTheta_Pos;
    TH1F*  fPhi_Pos;
    TH1F*  fVertAngle_Pos;
    TH1F*  fHorizAngle_Pos;
    TH2F*  fStartXvsStartY_Pos;
    TH2F*  fStartZvsStartX_Pos;
    TH2F*  fStartZvsStartY_Pos;
    TH2F*  fEndXvsEndY_Pos;
    TH2F*  fEndZvsEndX_Pos;
    TH2F*  fEndZvsEndY_Pos;
    TH1F*  fMinosX_Pos;
    TH1F*  fMinosY_Pos;
    TH1F*  fMinosZ_Pos;
    TH2F*  fMinosXY_Pos;
    TH1F*  fDirX_MINOS_start_Pos;
    TH1F*  fDirY_MINOS_start_Pos;
    TH1F*  fDirZ_MINOS_start_Pos;
    TH1F*  fCosX_MINOS_start_Pos;
    TH1F*  fCosY_MINOS_start_Pos;
    TH1F*  fCosZ_MINOS_start_Pos;
    TH1F*  fDiffDirX_Pos;
    TH1F*  fDiffDirY_Pos;
    TH1F*  fDiffDirZ_Pos;
    TH1F*  fDiffCosX_Pos;
    TH1F*  fDiffCosY_Pos;
    TH1F*  fDiffCosZ_Pos;
    TH1F*  fDiffX_Pos;
    TH1F*  fDiffY_Pos;
    TH1F*  fDiffR_Pos;

    TH2F*  fChannelVsHitAmplitude_Pos;
    //TH2F*  fChannelVsHitAmplitude_Corrected_Pos;

    TH1F*  fDirX_start_Neg;
    TH1F*  fDirY_start_Neg;
    TH1F*  fDirZ_start_Neg;
    TH1F*  fCosX_start_Neg;
    TH1F*  fCosY_start_Neg;
    TH1F*  fCosZ_start_Neg;
    TH1F*  fX_start_Neg;
    TH1F*  fY_start_Neg;
    TH1F*  fZ_start_Neg;
    TH1F*  fX_end_Neg;
    TH1F*  fY_end_Neg;
    TH1F*  fZ_end_Neg;
    TH1F*  fTrackLength_Neg;
    TH1F*  fTheta_Neg;
    TH1F*  fPhi_Neg;
    TH1F*  fVertAngle_Neg;
    TH1F*  fHorizAngle_Neg;
    TH2F*  fStartXvsStartY_Neg;
    TH2F*  fStartZvsStartX_Neg;
    TH2F*  fStartZvsStartY_Neg;
    TH2F*  fEndXvsEndY_Neg;
    TH2F*  fEndZvsEndX_Neg;
    TH2F*  fEndZvsEndY_Neg;
    TH1F*  fMinosX_Neg;
    TH1F*  fMinosY_Neg;
    TH1F*  fMinosZ_Neg;
    TH2F*  fMinosXY_Neg;
    TH1F*  fDirX_MINOS_start_Neg;
    TH1F*  fDirY_MINOS_start_Neg;
    TH1F*  fDirZ_MINOS_start_Neg;
    TH1F*  fCosX_MINOS_start_Neg;
    TH1F*  fCosY_MINOS_start_Neg;
    TH1F*  fCosZ_MINOS_start_Neg;
    TH1F*  fDiffDirX_Neg;
    TH1F*  fDiffDirY_Neg;
    TH1F*  fDiffDirZ_Neg;
    TH1F*  fDiffCosX_Neg;
    TH1F*  fDiffCosY_Neg;
    TH1F*  fDiffCosZ_Neg;
    TH1F*  fDiffX_Neg;
    TH1F*  fDiffY_Neg;
    TH1F*  fDiffR_Neg;

    TH2F*  fChannelVsHitAmplitude_Neg;
    //TH2F*  fChannelVsHitAmplitude_Corrected_Neg;

    TH1D* fMinosErange_Pos;
    TH1D* fMinosErange_Neg;
    TH1D* fMinosMom_Pos;
    TH1D* fMinosMom_Neg;
       
    TH1F* fMinosTrkChi2_Pos;
    TH1F* fMinosTrkChi2_Neg;
    TH2F*  fMinosTrkChi2vNPoints_Pos;
    TH2F*  fMinosTrkChi2vNPoints_Neg;

        
    std::string fTracks_label;
    std::string fTrackMatchModuleLabel;
    double fdBoundary; //distance from a boundary to be considered a track that "ends on a boundary"
        
  };
  
//-----------------------------------------------------------------------------
  MuonAna::MuonAna(fhicl::ParameterSet const& pset) :
    fTracks_label          (pset.get< std::string >("LArTracksModuleLabel")),
    fTrackMatchModuleLabel (pset.get< std::string >("TrackMatchModuleLabel")),
    fdBoundary             (pset.get< double >("dBoundary"))
  {
        
  }
    
//-----------------------------------------------------------------------------
  MuonAna::~MuonAna()
  {
  }

//-----------------------------------------------------------------------------
  void MuonAna::beginJob()
  {
    // get access to the TFile service  
    art::ServiceHandle<art::TFileService> tfs;

    fDirX_start = tfs->make<TH1F>("fDirX_start","Start Dir X", 180,0.0,180.0);
    fDirY_start = tfs->make<TH1F>("fDirY_start","Start Dir Y", 180,0.0,180.0);
    fDirZ_start = tfs->make<TH1F>("fDirZ_start","Start Dir Z", 180,0.0,180.0);
    fCosX_start = tfs->make<TH1F>("fCosX_start","Start Cos X", 100,-1.0,1.0);
    fCosY_start = tfs->make<TH1F>("fCosY_start","Start Cos Y", 100,-1.0,1.0);
    fCosZ_start = tfs->make<TH1F>("fCosZ_start","Start Cos Z", 100,-1.0,1.0);
    fX_start = tfs->make<TH1F>("fX_start","Start X", 220,-5.0,50.0);
    fY_start = tfs->make<TH1F>("fY_start","Start Y", 200,-25.0,25.0);
    fZ_start = tfs->make<TH1F>("fZ_start","Start Z", 400,-5.0,95.0);
    fX_end = tfs->make<TH1F>("fX_end","End X", 220,-5.0,50.0);
    fY_end = tfs->make<TH1F>("fY_end","End Y", 200,-25.0,25.0);
    fZ_end = tfs->make<TH1F>("fZ_end","End Z", 400,-5.0,95.0);
    fTrackLength = tfs->make<TH1F>("fTrackLength","TrackLength", 600,0.0,150.0);
    fTheta = tfs->make<TH1F>("fTheta","Theta", 100,0.0,TMath::Pi());
    fPhi = tfs->make<TH1F>("fPhi","Phi", 200,-1.0*TMath::Pi(),TMath::Pi());
    fVertAngle = tfs->make<TH1F>("fVertAngle","VertAngle", 180,-90.0,90.0);
    fHorizAngle = tfs->make<TH1F>("fHorizAngle","HorizAngle", 180,-90.0,90.0);
    fStartXvsStartY = tfs->make<TH2F>("fStartXvsStartY","Start X vs. Y", 220,-5.0,50.0,200,-25.0,25.0);
    fStartZvsStartX = tfs->make<TH2F>("fStartZvsStartX","Start Z vs. X", 400,-5.0,95.0,220,-5.0,50.0);
    fStartZvsStartY = tfs->make<TH2F>("fStartZvsStartY","Start Z vs. Y", 400,-5.0,95.0,200,-25.0,25.0);
    fEndXvsEndY = tfs->make<TH2F>("fEndXvsEndY","End X vs. Y", 220,-5.0,50.0,200,-25.0,25.0);
    fEndZvsEndX = tfs->make<TH2F>("fEndZvsEndX","End Z vs. X", 400,-5.0,95.0,220,-5.0,50.0);
    fEndZvsEndY = tfs->make<TH2F>("fEndZvsEndY","End Z vs. Y", 400,-5.0,95.0,200,-25.0,25.0);
    fChannelVsHitAmplitude = tfs->make<TH2F>("fChannelVsHitAmplitude","Channel vs. Hit Amplitude",480,0.0,480.0,1000,0.0,200.0);
    //fChannelVsHitAmplitude_Corrected = tfs->make<TH2F>("fChannelVsHitAmplitude_Corrected","Channel vs. Hit Amplitude, Lifetime",480,0.0,480.0,1000,0.0,200.0);

    //TG Muons matched to Pos./Neg. MINOS tracks
    fDirX_start_Pos = tfs->make<TH1F>("fDirX_start_Pos","Start Dir X Pos", 180,0.0,180.0);
    fDirY_start_Pos = tfs->make<TH1F>("fDirY_start_Pos","Start Dir Y Pos", 180,0.0,180.0);
    fDirZ_start_Pos = tfs->make<TH1F>("fDirZ_start_Pos","Start Dir Z Pos", 180,0.0,180.0);
    fCosX_start_Pos = tfs->make<TH1F>("fCosX_start_Pos","Start Cos X Pos", 100,-1.0,1.0);
    fCosY_start_Pos = tfs->make<TH1F>("fCosY_start_Pos","Start Cos Y Pos", 100,-1.0,1.0);
    fCosZ_start_Pos = tfs->make<TH1F>("fCosZ_start_Pos","Start Cos Z Pos", 100,-1.0,1.0);
    fX_start_Pos = tfs->make<TH1F>("fX_start_Pos","Start X Pos", 220,-5.0,50.0);
    fY_start_Pos = tfs->make<TH1F>("fY_start_Pos","Start Y Pos", 200,-25.0,25.0);
    fZ_start_Pos = tfs->make<TH1F>("fZ_start_Pos","Start Z Pos", 400,-5.0,95.0);
    fX_end_Pos = tfs->make<TH1F>("fX_end_Pos","End X Pos", 220,-5.0,50.0);
    fY_end_Pos = tfs->make<TH1F>("fY_end_Pos","End Y Pos", 200,-25.0,25.0);
    fZ_end_Pos = tfs->make<TH1F>("fZ_end_Pos","End Z Pos", 400,-5.0,95.0);
    fTrackLength_Pos = tfs->make<TH1F>("fTrackLength_Pos","TrackLength Pos", 600,0.0,150.0);
    fTheta_Pos = tfs->make<TH1F>("fTheta_Pos","Theta Pos", 100,0.0,TMath::Pi());
    fPhi_Pos = tfs->make<TH1F>("fPhi_Pos","Phi Pos", 200,-1.0*TMath::Pi(),TMath::Pi());
    fVertAngle_Pos = tfs->make<TH1F>("fVertAngle_Pos","VertAngle Pos", 180,-90.0,90.0);
    fHorizAngle_Pos = tfs->make<TH1F>("fHorizAngle_Pos","HorizAngle Pos", 180,-90.0,90.0);
    fStartXvsStartY_Pos = tfs->make<TH2F>("fStartXvsStartY_Pos","Start X vs. Y Pos", 220,-5.0,50.0,200,-25.0,25.0);
    fStartZvsStartX_Pos = tfs->make<TH2F>("fStartZvsStartX_Pos","Start Z vs. X Pos", 400,-5.0,95.0,220,-5.0,50.0);
    fStartZvsStartY_Pos = tfs->make<TH2F>("fStartZvsStartY_Pos","Start Z vs. Y Pos", 400,-5.0,95.0,200,-25.0,25.0);
    fEndXvsEndY_Pos = tfs->make<TH2F>("fEndXvsEndY_Pos","End X vs. Y Pos", 220,-5.0,50.0,200,-25.0,25.0);
    fEndZvsEndX_Pos = tfs->make<TH2F>("fEndZvsEndX_Pos","End Z vs. X Pos", 400,-5.0,95.0,220,-5.0,50.0);
    fEndZvsEndY_Pos = tfs->make<TH2F>("fEndZvsEndY_Pos","End Z vs. Y Pos", 400,-5.0,95.0,200,-25.0,25.0);
    fMinosX_Pos = tfs->make<TH1F>("fMinosX_Pos","Minos Vtx X Pos",600,-300.0,300.0);
    fMinosY_Pos = tfs->make<TH1F>("fMinosY_Pos","Minos Vtx Y Pos",300,-150.0,150.0);
    fMinosZ_Pos = tfs->make<TH1F>("fMinosZ_Pos","Minos Vtx Z Pos",150,-50.0,100.0);
    fMinosXY_Pos = tfs->make<TH2F>("fMinosXY_Pos","Minos Vtx X vs. Y Pos",600,-300.0,300.0,300,-150.0,150.0);
    fDirX_MINOS_start_Pos = tfs->make<TH1F>("fDirX_MINOS_start_Pos","MINOS Start Dir X Pos", 180,0.0,180.0);
    fDirY_MINOS_start_Pos = tfs->make<TH1F>("fDirY_MINOS_start_Pos","MINOS Start Dir Y Pos", 180,0.0,180.0);
    fDirZ_MINOS_start_Pos = tfs->make<TH1F>("fDirZ_MINOS_start_Pos","MINOS Start Dir Z Pos", 180,0.0,180.0);
    fCosX_MINOS_start_Pos = tfs->make<TH1F>("fCosX_MINOS_start_Pos","MINOS Start Cos X Pos", 100,-1.0,1.0);
    fCosY_MINOS_start_Pos = tfs->make<TH1F>("fCosY_MINOS_start_Pos","MINOS Start Cos Y Pos", 100,-1.0,1.0);
    fCosZ_MINOS_start_Pos = tfs->make<TH1F>("fCosZ_MINOS_start_Pos","MINOS Start Cos Z Pos", 100,-1.0,1.0);
    fDiffDirX_Pos = tfs->make<TH1F>("fDiffDirX_Pos","Diff Dir X Pos", 720,-180.0,180.0);
    fDiffDirY_Pos = tfs->make<TH1F>("fDiffDirY_Pos","Diff Dir Y Pos", 720,-180.0,180.0);
    fDiffDirZ_Pos = tfs->make<TH1F>("fDiffDirZ_Pos","Diff Dir Z Pos", 720,-180.0,180.0);
    fDiffCosX_Pos = tfs->make<TH1F>("fDiffCosX_Pos","Diff Cos X Pos", 800,-2.0,2.0);
    fDiffCosY_Pos = tfs->make<TH1F>("fDiffCosY_Pos","Diff Cos Y Pos", 800,-2.0,2.0);
    fDiffCosZ_Pos = tfs->make<TH1F>("fDiffCosZ_Pos","Diff Cos Z Pos", 800,-2.0,2.0);
    fDiffX_Pos = tfs->make<TH1F>("fDiffX_Pos","Diff X Pos", 200, -100.0, 100.0);
    fDiffY_Pos = tfs->make<TH1F>("fDiffY_Pos","Diff Y Pos", 200, -100.0, 100.0);
    fDiffR_Pos = tfs->make<TH1F>("fDiffR_Pos","Diff R Pos", 200, 0.0, 200.0);
    fChannelVsHitAmplitude_Pos = tfs->make<TH2F>("fChannelVsHitAmplitude_Pos","Channel vs. Hit Amplitude Pos",480,0.0,480.0,1000,0.0,200.0);
    //fChannelVsHitAmplitude_Corrected_Pos = tfs->make<TH2F>("fChannelVsHitAmplitude_Corrected_Pos","Channel vs. Hit Amplitude, Lifetime, Pos",480,0.0,480.0,1000,0.0,200.0);

    fDirX_start_Neg = tfs->make<TH1F>("fDirX_start_Neg","Start Dir X Neg", 180,0.0,180.0);
    fDirY_start_Neg = tfs->make<TH1F>("fDirY_start_Neg","Start Dir Y Neg", 180,0.0,180.0);
    fDirZ_start_Neg = tfs->make<TH1F>("fDirZ_start_Neg","Start Dir Z Neg", 180,0.0,180.0);
    fCosX_start_Neg = tfs->make<TH1F>("fCosX_start_Neg","Start Cos X Neg", 100,-1.0,1.0);
    fCosY_start_Neg = tfs->make<TH1F>("fCosY_start_Neg","Start Cos Y Neg", 100,-1.0,1.0);
    fCosZ_start_Neg = tfs->make<TH1F>("fCosZ_start_Neg","Start Cos Z Neg", 100,-1.0,1.0);
    fX_start_Neg = tfs->make<TH1F>("fX_start_Neg","Start X Neg", 220,-5.0,50.0);
    fY_start_Neg = tfs->make<TH1F>("fY_start_Neg","Start Y Neg", 200,-25.0,25.0);
    fZ_start_Neg = tfs->make<TH1F>("fZ_start_Neg","Start Z Neg", 400,-5.0,95.0);
    fX_end_Neg = tfs->make<TH1F>("fX_end_Neg","End X Neg", 220,-5.0,50.0);
    fY_end_Neg = tfs->make<TH1F>("fY_end_Neg","End Y Neg", 200,-25.0,25.0);
    fZ_end_Neg = tfs->make<TH1F>("fZ_end_Neg","End Z Neg", 400,-5.0,95.0);
    fTrackLength_Neg = tfs->make<TH1F>("fTrackLength_Neg","TrackLength Neg", 600,0.0,150.0);
    fTheta_Neg = tfs->make<TH1F>("fTheta_Neg","Theta Neg", 100,0.0,TMath::Pi());
    fPhi_Neg = tfs->make<TH1F>("fPhi_Neg","Phi Neg", 200,-1.0*TMath::Pi(),TMath::Pi());
    fVertAngle_Neg = tfs->make<TH1F>("fVertAngle_Neg","VertAngle Neg", 180,-90.0,90.0);
    fHorizAngle_Neg = tfs->make<TH1F>("fHorizAngle_Neg","HorizAngle Neg", 180,-90.0,90.0);
    fStartXvsStartY_Neg = tfs->make<TH2F>("fStartXvsStartY_Neg","Start X vs. Y Neg", 220,-5.0,50.0,200,-25.0,25.0);
    fStartZvsStartX_Neg = tfs->make<TH2F>("fStartZvsStartX_Neg","Start Z vs. X Neg", 400,-5.0,95.0,220,-5.0,50.0);
    fStartZvsStartY_Neg = tfs->make<TH2F>("fStartZvsStartY_Neg","Start Z vs. Y Neg", 400,-5.0,95.0,200,-25.0,25.0);
    fEndXvsEndY_Neg = tfs->make<TH2F>("fEndXvsEndY_Neg","End X vs. Y Neg", 220,-5.0,50.0,200,-25.0,25.0);
    fEndZvsEndX_Neg = tfs->make<TH2F>("fEndZvsEndX_Neg","End Z vs. X Neg", 400,-5.0,95.0,220,-5.0,50.0);
    fEndZvsEndY_Neg = tfs->make<TH2F>("fEndZvsEndY_Neg","End Z vs. Y Neg", 400,-5.0,95.0,200,-25.0,25.0);
    fMinosX_Neg = tfs->make<TH1F>("fMinosX_Neg","Minos Vtx X Neg",600,-300.0,300.0);
    fMinosY_Neg = tfs->make<TH1F>("fMinosY_Neg","Minos Vtx Y Neg",300,-150.0,150.0);
    fMinosZ_Neg = tfs->make<TH1F>("fMinosZ_Neg","Minos Vtx Z Neg",150,-50.0,100.0);
    fMinosXY_Neg = tfs->make<TH2F>("fMinosXY_Neg","Minos Vtx X vs. Y Neg",600,-300.0,300.0,300,-150.0,150.0);
    fDirX_MINOS_start_Neg = tfs->make<TH1F>("fDirX_MINOS_start_Neg","MINOS Start Dir X Neg", 180,0.0,180.0);
    fDirY_MINOS_start_Neg = tfs->make<TH1F>("fDirY_MINOS_start_Neg","MINOS Start Dir Y Neg", 180,0.0,180.0);
    fDirZ_MINOS_start_Neg = tfs->make<TH1F>("fDirZ_MINOS_start_Neg","MINOS Start Dir Z Neg", 180,0.0,180.0);
    fCosX_MINOS_start_Neg = tfs->make<TH1F>("fCosX_MINOS_start_Neg","MINOS Start Cos X Neg", 100,-1.0,1.0);
    fCosY_MINOS_start_Neg = tfs->make<TH1F>("fCosY_MINOS_start_Neg","MINOS Start Cos Y Neg", 100,-1.0,1.0);
    fCosZ_MINOS_start_Neg = tfs->make<TH1F>("fCosZ_MINOS_start_Neg","MINOS Start Cos Z Neg", 100,-1.0,1.0);
    fDiffDirX_Neg = tfs->make<TH1F>("fDiffDirX_Neg","Diff Dir X Neg", 720,-180.0,180.0);
    fDiffDirY_Neg = tfs->make<TH1F>("fDiffDirY_Neg","Diff Dir Y Neg", 720,-180.0,180.0);
    fDiffDirZ_Neg = tfs->make<TH1F>("fDiffDirZ_Neg","Diff Dir Z Neg", 720,-180.0,180.0);
    fDiffCosX_Neg = tfs->make<TH1F>("fDiffCosX_Neg","Diff Cos X Neg", 800,-2.0,2.0);
    fDiffCosY_Neg = tfs->make<TH1F>("fDiffCosY_Neg","Diff Cos Y Neg", 800,-2.0,2.0);
    fDiffCosZ_Neg = tfs->make<TH1F>("fDiffCosZ_Neg","Diff Cos Z Neg", 800,-2.0,2.0);
    fDiffX_Neg = tfs->make<TH1F>("fDiffX_Neg","Diff X Neg", 200, -100.0, 100.0);
    fDiffY_Neg = tfs->make<TH1F>("fDiffY_Neg","Diff Y Neg", 200, -100.0, 100.0);
    fDiffR_Neg = tfs->make<TH1F>("fDiffR_Neg","Diff R Neg", 200, 0.0, 200.0);
    fChannelVsHitAmplitude_Neg = tfs->make<TH2F>("fChannelVsHitAmplitude_Neg","Channel vs. Hit Amplitude Neg",480,0.0,480.0,1000,0.0,200.0);
    //fChannelVsHitAmplitude_Corrected_Neg = tfs->make<TH2F>("fChannelVsHitAmplitude_Corrected_Neg","Channel vs. Hit Amplitude, Lifetime, Neg",480,0.0,480.0,1000,0.0,200.0);

    fMinosErange_Pos = tfs->make<TH1D>("fMinosErange_Pos","MINOS + Charge Tracks: Erange",5000,0.0,50.0);
    fMinosErange_Neg = tfs->make<TH1D>("fMinosErange_Neg","MINOS - Charge Tracks: Erange",5000,0.0,50.0);
    fMinosMom_Pos = tfs->make<TH1D>("fMinosMom_Pos","MINOS + Charge Tracks: Momentum",5000,0.0,50.0);
    fMinosMom_Neg = tfs->make<TH1D>("fMinosMom_Neg","MINOS - Charge Tracks: Momentum",5000,0.0,50.0);

    fMinosTrkChi2_Pos = tfs->make<TH1F>("fMinosTrkChi2_Pos","MINOS TrkChi2 Pos",100,0.0,2000.0);
    fMinosTrkChi2_Neg = tfs->make<TH1F>("fMinosTrkChi2_Neg","MINOS TrkChi2 Neg",100,0.0,2000.0);

    fMinosTrkChi2vNPoints_Pos = tfs->make<TH2F>("fMinosTrkChi2vNPoints_Pos","MINOS TrkChi2 Pos vs. Points",100,0.0,2000.0, 40, 0.0, 200.0);
    fMinosTrkChi2vNPoints_Neg = tfs->make<TH2F>("fMinosTrkChi2vNPoints_Neg","MINOS TrkChi2 Neg vs. Points",100,0.0,2000.0, 40, 0.0, 200.0);
        
  }
    
//-----------------------------------------------------------------------------
  void MuonAna::analyze(const art::Event& evt) 
  {
    art::Handle< std::vector<recob::Track> > LarTrackHandle;
    evt.getByLabel(fTracks_label,LarTrackHandle);
        
    double trackCosStart[3]={0.,0.,0.};
    double trackCosEnd[3]={0.,0.,0.};
    std::vector<double> larStart;
    std::vector<double> larEnd;
        
    //find matched MINOS information for each track
    art::FindOne<t962::MINOS> fomatch(LarTrackHandle, evt, fTrackMatchModuleLabel);

    //find Hit information for each track
    art::FindManyP<recob::Hit> fhit(LarTrackHandle, evt, fTracks_label);

      
    for(unsigned int i=0; i<LarTrackHandle->size();++i){
                
      art::Ptr<recob::Track> lartrack(LarTrackHandle,i);
              
      bool startsonboundary = BeginsOnBoundary(lartrack);
      bool endsonboundary   = EndsOnBoundary(lartrack);

      if(startsonboundary && endsonboundary){//throughgoing track...make plots
        lartrack->Extent(larStart,larEnd);
        lartrack->Direction(trackCosStart,trackCosEnd);

        double tracklength = sqrt((larStart[0]-larEnd[0])*(larStart[0]-larEnd[0]) +
                                  (larStart[1]-larEnd[1])*(larStart[1]-larEnd[1]) +
                                  (larStart[2]-larEnd[2])*(larStart[2]-larEnd[2]));

        fDirX_start->Fill(TMath::ACos(trackCosStart[0])*180.0/TMath::Pi());
        fDirY_start->Fill(TMath::ACos(trackCosStart[1])*180.0/TMath::Pi());
        fDirZ_start->Fill(TMath::ACos(trackCosStart[2])*180.0/TMath::Pi());
        fCosX_start->Fill(trackCosStart[0]);
        fCosY_start->Fill(trackCosStart[1]);
        fCosZ_start->Fill(trackCosStart[2]);
        fX_start->Fill(larStart[0]);
        fY_start->Fill(larStart[1]);
        fZ_start->Fill(larStart[2]);
        fX_end->Fill(larEnd[0]);
        fY_end->Fill(larEnd[1]);
        fZ_end->Fill(larEnd[2]);
        fTrackLength->Fill(tracklength);
        fTheta->Fill(lartrack->Theta());
        fPhi->Fill(lartrack->Phi());
        fVertAngle->Fill(TMath::ATan(trackCosStart[1]/trackCosStart[2])*180.0/TMath::Pi());
        fHorizAngle->Fill(TMath::ATan(trackCosStart[0]/trackCosStart[2])*180.0/TMath::Pi());
        fStartXvsStartY->Fill(larStart[0],larStart[1]);
        fStartZvsStartX->Fill(larStart[2],larStart[0]);
        fStartZvsStartY->Fill(larStart[2],larStart[1]);
        fEndXvsEndY->Fill(larEnd[0],larEnd[1]);
        fEndZvsEndX->Fill(larEnd[2],larEnd[0]);
        fEndZvsEndY->Fill(larEnd[2],larEnd[1]);

        if(fhit.at(i).size()>0){
          std::vector< art::Ptr<recob::Hit> > trackhits = fhit.at(i);
          if(trackhits.size()>100){
            for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = trackhits.begin();
                ihit != trackhits.end(); ++ihit){
              const recob::Hit& hit = **ihit;
              fChannelVsHitAmplitude->Fill(hit.Channel()+0.1,hit.Charge(true));
              //fChannelVsHitAmplitude->Fill(hit.Channel(),hit.Charge(true));
            }
          }
        }



        if(!fomatch.at(i).isValid()) continue;//No matching MINOS track

        double xdiff, ydiff, rdiff;
        Compare(lartrack, fomatch.at(i).ref(), xdiff, ydiff, rdiff);

        if(fomatch.at(i).ref().fcharge==1){
          fDirX_start_Pos->Fill(TMath::ACos(trackCosStart[0])*180.0/TMath::Pi());
          fDirY_start_Pos->Fill(TMath::ACos(trackCosStart[1])*180.0/TMath::Pi());
          fDirZ_start_Pos->Fill(TMath::ACos(trackCosStart[2])*180.0/TMath::Pi());
          fCosX_start_Pos->Fill(trackCosStart[0]);
          fCosY_start_Pos->Fill(trackCosStart[1]);
          fCosZ_start_Pos->Fill(trackCosStart[2]);
          fX_start_Pos->Fill(larStart[0]);
          fY_start_Pos->Fill(larStart[1]);
          fZ_start_Pos->Fill(larStart[2]);
          fX_end_Pos->Fill(larEnd[0]);
          fY_end_Pos->Fill(larEnd[1]);
          fZ_end_Pos->Fill(larEnd[2]);
          fTrackLength_Pos->Fill(tracklength);
          fTheta_Pos->Fill(lartrack->Theta());
          fPhi_Pos->Fill(lartrack->Phi());
          fVertAngle_Pos->Fill(TMath::ATan(trackCosStart[1]/trackCosStart[2])*180.0/TMath::Pi());
          fHorizAngle_Pos->Fill(TMath::ATan(trackCosStart[0]/trackCosStart[2])*180.0/TMath::Pi());
          fStartXvsStartY_Pos->Fill(larStart[0],larStart[1]);
          fStartZvsStartX_Pos->Fill(larStart[2],larStart[0]);
          fStartZvsStartY_Pos->Fill(larStart[2],larStart[1]);
          fEndXvsEndY_Pos->Fill(larEnd[0],larEnd[1]);
          fEndZvsEndX_Pos->Fill(larEnd[2],larEnd[0]);
          fEndZvsEndY_Pos->Fill(larEnd[2],larEnd[1]);
          fMinosX_Pos->Fill(100.0*fomatch.at(i).ref().ftrkVtxX);
          fMinosY_Pos->Fill(100.0*fomatch.at(i).ref().ftrkVtxY);
          fMinosZ_Pos->Fill(100.0*fomatch.at(i).ref().ftrkVtxZ);
          fMinosXY_Pos->Fill(100.0*fomatch.at(i).ref().ftrkVtxX,100.0*fomatch.at(i).ref().ftrkVtxY);
          fDirX_MINOS_start_Pos->Fill(TMath::ACos(fomatch.at(i).ref().ftrkdcosx)*180.0/TMath::Pi());
          fDirY_MINOS_start_Pos->Fill(TMath::ACos(fomatch.at(i).ref().ftrkdcosy)*180.0/TMath::Pi());
          fDirZ_MINOS_start_Pos->Fill(TMath::ACos(fomatch.at(i).ref().ftrkdcosz)*180.0/TMath::Pi());
          fCosX_MINOS_start_Pos->Fill(fomatch.at(i).ref().ftrkdcosx);
          fCosY_MINOS_start_Pos->Fill(fomatch.at(i).ref().ftrkdcosy);
          fCosZ_MINOS_start_Pos->Fill(fomatch.at(i).ref().ftrkdcosz);
          fDiffDirX_Pos->Fill((TMath::ACos(trackCosStart[0]) - TMath::ACos(fomatch.at(i).ref().ftrkdcosx))*180.0/TMath::Pi());
          fDiffDirY_Pos->Fill((TMath::ACos(trackCosStart[1]) - TMath::ACos(fomatch.at(i).ref().ftrkdcosy))*180.0/TMath::Pi());
          fDiffDirZ_Pos->Fill((TMath::ACos(trackCosStart[2]) - TMath::ACos(fomatch.at(i).ref().ftrkdcosz))*180.0/TMath::Pi());
          fDiffCosX_Pos->Fill(  trackCosStart[0] - fomatch.at(i).ref().ftrkdcosx);
          fDiffCosY_Pos->Fill(  trackCosStart[1] - fomatch.at(i).ref().ftrkdcosy);
          fDiffCosZ_Pos->Fill(  trackCosStart[2] - fomatch.at(i).ref().ftrkdcosz);
          fDiffX_Pos->Fill(xdiff);
          fDiffY_Pos->Fill(ydiff);
          fDiffR_Pos->Fill(rdiff);
               
          if(fhit.at(i).size()>0){
            std::vector< art::Ptr<recob::Hit> > trackhits = fhit.at(i);
            if(trackhits.size()>100){
              for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = trackhits.begin();
                  ihit != trackhits.end(); ++ihit){
                const recob::Hit& hit = **ihit;
                fChannelVsHitAmplitude_Pos->Fill(hit.Channel()+0.1,hit.Charge(true));
                //fChannelVsHitAmplitude_Pos->Fill(hit.Channel(),hit.Charge(true));
              }
            }
          }
               


          fMinosTrkChi2_Pos->Fill(fomatch.at(i).ref().ftrkChi2); 
          fMinosTrkChi2vNPoints_Pos->Fill(fomatch.at(i).ref().ftrkChi2,fomatch.at(i).ref().fntrkstp);
                           
                        
          if(fomatch.at(i).ref().ftrkcontained){
            fMinosMom_Pos->Fill(fomatch.at(i).ref().ftrkErange);
            fMinosErange_Pos->Fill(fomatch.at(i).ref().ftrkErange);
          }
          else fMinosMom_Pos->Fill(fomatch.at(i).ref().ftrkmom);
        }
        if(fomatch.at(i).ref().fcharge==-1){
          fDirX_start_Neg->Fill(TMath::ACos(trackCosStart[0])*180.0/TMath::Pi());
          fDirY_start_Neg->Fill(TMath::ACos(trackCosStart[1])*180.0/TMath::Pi());
          fDirZ_start_Neg->Fill(TMath::ACos(trackCosStart[2])*180.0/TMath::Pi());
          fCosX_start_Neg->Fill(trackCosStart[0]);
          fCosY_start_Neg->Fill(trackCosStart[1]);
          fCosZ_start_Neg->Fill(trackCosStart[2]);
          fX_start_Neg->Fill(larStart[0]);
          fY_start_Neg->Fill(larStart[1]);
          fZ_start_Neg->Fill(larStart[2]);
          fX_end_Neg->Fill(larEnd[0]);
          fY_end_Neg->Fill(larEnd[1]);
          fZ_end_Neg->Fill(larEnd[2]);
          fTrackLength_Neg->Fill(tracklength);
          fTheta_Neg->Fill(lartrack->Theta());
          fPhi_Neg->Fill(lartrack->Phi());
          fVertAngle_Neg->Fill(TMath::ATan(trackCosStart[1]/trackCosStart[2])*180.0/TMath::Pi());
          fHorizAngle_Neg->Fill(TMath::ATan(trackCosStart[0]/trackCosStart[2])*180.0/TMath::Pi());
          fStartXvsStartY_Neg->Fill(larStart[0],larStart[1]);
          fStartZvsStartX_Neg->Fill(larStart[2],larStart[0]);
          fStartZvsStartY_Neg->Fill(larStart[2],larStart[1]);
          fEndXvsEndY_Neg->Fill(larEnd[0],larEnd[1]);
          fEndZvsEndX_Neg->Fill(larEnd[2],larEnd[0]);
          fEndZvsEndY_Neg->Fill(larEnd[2],larEnd[1]);
          fMinosX_Neg->Fill(100.0*fomatch.at(i).ref().ftrkVtxX);
          fMinosY_Neg->Fill(100.0*fomatch.at(i).ref().ftrkVtxY);
          fMinosZ_Neg->Fill(100.0*fomatch.at(i).ref().ftrkVtxZ);
          fMinosXY_Neg->Fill(100.0*fomatch.at(i).ref().ftrkVtxX,100.0*fomatch.at(i).ref().ftrkVtxY);
          fDirX_MINOS_start_Neg->Fill(TMath::ACos(fomatch.at(i).ref().ftrkdcosx)*180.0/TMath::Pi());
          fDirY_MINOS_start_Neg->Fill(TMath::ACos(fomatch.at(i).ref().ftrkdcosy)*180.0/TMath::Pi());
          fDirZ_MINOS_start_Neg->Fill(TMath::ACos(fomatch.at(i).ref().ftrkdcosz)*180.0/TMath::Pi());
          fCosX_MINOS_start_Neg->Fill(fomatch.at(i).ref().ftrkdcosx);
          fCosY_MINOS_start_Neg->Fill(fomatch.at(i).ref().ftrkdcosy);
          fCosZ_MINOS_start_Neg->Fill(fomatch.at(i).ref().ftrkdcosz);
          fDiffDirX_Neg->Fill((TMath::ACos(trackCosStart[0]) - TMath::ACos(fomatch.at(i).ref().ftrkdcosx))*180.0/TMath::Pi());
          fDiffDirY_Neg->Fill((TMath::ACos(trackCosStart[1]) - TMath::ACos(fomatch.at(i).ref().ftrkdcosy))*180.0/TMath::Pi());
          fDiffDirZ_Neg->Fill((TMath::ACos(trackCosStart[2]) - TMath::ACos(fomatch.at(i).ref().ftrkdcosz))*180.0/TMath::Pi());
          fDiffCosX_Neg->Fill(  trackCosStart[0] - fomatch.at(i).ref().ftrkdcosx);
          fDiffCosY_Neg->Fill(  trackCosStart[1] - fomatch.at(i).ref().ftrkdcosy);
          fDiffCosZ_Neg->Fill(  trackCosStart[2] - fomatch.at(i).ref().ftrkdcosz);
          fDiffX_Neg->Fill(xdiff);
          fDiffY_Neg->Fill(ydiff);
          fDiffR_Neg->Fill(rdiff);

          if(fhit.at(i).size()>0){
            std::vector< art::Ptr<recob::Hit> > trackhits = fhit.at(i);
            if(trackhits.size()>100){
              for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = trackhits.begin();
                  ihit != trackhits.end(); ++ihit){
                const recob::Hit& hit = **ihit;
                fChannelVsHitAmplitude_Neg->Fill(hit.Channel()+0.1,hit.Charge(true));
                //fChannelVsHitAmplitude_Neg->Fill(hit.Channel(),hit.Charge(true));
              }
            }
          }
               

          fMinosTrkChi2_Neg->Fill(fomatch.at(i).ref().ftrkChi2); 
          fMinosTrkChi2vNPoints_Neg->Fill(fomatch.at(i).ref().ftrkChi2,fomatch.at(i).ref().fntrkstp);
                           
          if(fomatch.at(i).ref().ftrkcontained){
            fMinosMom_Neg->Fill(fomatch.at(i).ref().ftrkErange);
            fMinosErange_Neg->Fill(fomatch.at(i).ref().ftrkErange);
          }
          else fMinosMom_Neg->Fill(fomatch.at(i).ref().ftrkmom);
        }
         
      }//through-going 
                
    }//loop over Tracks

  }//analyze

//--------------------------------------------------
  bool MuonAna::BeginsOnBoundary(art::Ptr<recob::Track> lar_track)
  {
    std::vector<double> larStart, larEnd;
    lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)
    if(fabs(larStart[0])<fdBoundary
       || fabs(47.-larStart[0])<fdBoundary 
       || fabs(larStart[1]+20.)<fdBoundary
       || fabs(20.-larStart[1])<fdBoundary
       || fabs(larStart[2])<fdBoundary
       || fabs(90.-larStart[2])<fdBoundary )   
      return true;  
    else return false;
  }
    
//--------------------------------------------------
  bool MuonAna::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
  {
    std::vector<double> larStart, larEnd;
    lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)
    if(fabs(larEnd[0])<fdBoundary
       || fabs(47.-larEnd[0])<fdBoundary 
       || fabs(larEnd[1]+20.)<fdBoundary
       || fabs(20.-larEnd[1])<fdBoundary
       || fabs(larEnd[2])<fdBoundary
       || fabs(90.-larEnd[2])<fdBoundary )   
      return true;  
    else return false;
  }

//--------------------------------------------------
  void MuonAna::Compare(art::Ptr<recob::Track> lar_track, t962::MINOS minos_track,
                        double &dx, double &dy, double &rdiff)
  {
    double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane 
    //(this minus number is the one we measured with Mitch)

    double x_offset=117.4; // previously 116.9;
    double y_offset=19.3; // previously  20.28;
    std::vector<double> larStart, larEnd;
    lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)

    double lardirectionStart[3];
    double lardirectionEnd[3];
    lar_track->Direction(lardirectionStart,lardirectionEnd);


    double dz = D - larEnd[2]+(100.0 * minos_track.ftrkVtxZ);//z-difference between end of T962 track and
    //begin of MINOS track...in centimeters

    double l = dz/(lardirectionEnd[2]);//3-d distance between end of T962 track and begin of MINOS track

    double x_pred = l*lardirectionEnd[0]+larEnd[0];//predicted x-pos. of T962 track at z-position equal to
    //start of MINOS track
    double y_pred = l*lardirectionEnd[1]+larEnd[1];//predicted y-pos. of T962 track at z-position equal to
    //start of MINOS track

    dx = 100.0*minos_track.ftrkVtxX - x_offset - x_pred;
    dy = 100.0*minos_track.ftrkVtxY + y_offset - y_pred;

    rdiff = sqrt(dx*dx + dy*dy);

    return;

  }

  //Required
  DEFINE_ART_MODULE(MuonAna);

} // end of namespace

#endif // MUONANA_H
