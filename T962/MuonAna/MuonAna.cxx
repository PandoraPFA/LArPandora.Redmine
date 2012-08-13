////////////////////////////////////////////////////////////////////////
//
// MuonAna class
//
// msoderbe@syr.edu
//
// Make plots for through-going muons.
////////////////////////////////////////////////////////////////////////


#include <iostream>


// Framework includes
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

#include "T962/MuonAna/MuonAna.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"

#include "Filters/ChannelFilter.h"
#include "T962/T962_Objects/MINOS.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

namespace muons {
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
      fDiffDirX_Pos = tfs->make<TH1F>("fDiffDirX_Pos","Diff Dir X Pos", 360,-180.0,180.0);
      fDiffDirY_Pos = tfs->make<TH1F>("fDiffDirY_Pos","Diff Dir Y Pos", 360,-180.0,180.0);
      fDiffDirZ_Pos = tfs->make<TH1F>("fDiffDirZ_Pos","Diff Dir Z Pos", 360,-180.0,180.0);
      fDiffCosX_Pos = tfs->make<TH1F>("fDiffCosX_Pos","Diff Cos X Pos", 400,-2.0,2.0);
      fDiffCosY_Pos = tfs->make<TH1F>("fDiffCosY_Pos","Diff Cos Y Pos", 400,-2.0,2.0);
      fDiffCosZ_Pos = tfs->make<TH1F>("fDiffCosZ_Pos","Diff Cos Z Pos", 400,-2.0,2.0);
      
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
      fDiffDirX_Neg = tfs->make<TH1F>("fDiffDirX_Neg","Diff Dir X Neg", 360,-180.0,180.0);
      fDiffDirY_Neg = tfs->make<TH1F>("fDiffDirY_Neg","Diff Dir Y Neg", 360,-180.0,180.0);
      fDiffDirZ_Neg = tfs->make<TH1F>("fDiffDirZ_Neg","Diff Dir Z Neg", 360,-180.0,180.0);
      fDiffCosX_Neg = tfs->make<TH1F>("fDiffCosX_Neg","Diff Cos X Neg", 400,-2.0,2.0);
      fDiffCosY_Neg = tfs->make<TH1F>("fDiffCosY_Neg","Diff Cos Y Neg", 400,-2.0,2.0);
      fDiffCosZ_Neg = tfs->make<TH1F>("fDiffCosZ_Neg","Diff Cos Z Neg", 400,-2.0,2.0);

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

            if(!fomatch.at(i).isValid()) continue;//No matching MINOS track

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


}//namespace
