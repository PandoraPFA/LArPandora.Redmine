////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to analyze Muons
///
/// \author  msoderbe@syr.edu
////////////////////////////////////////////////////////////////////////

#ifndef MUONANA_H
#define MUONANA_H

#include "art/Persistency/Common/Ptr.h"
#include "RecoBase/Track.h"

#include "Utilities/LArProperties.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>

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
        
    private:
       
        //plots for TG Muons
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

       //plots for TG Muons matched to a pos/neg MINOS track
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

        TH1D* fMinosErange_Pos;
        TH1D* fMinosErange_Neg;
        TH1D* fMinosMom_Pos;
        TH1D* fMinosMom_Neg;
       
       TH1F* fMinosTrkChi2_Pos;
       TH1F* fMinosTrkChi2_Neg;
       TH2F*  fMinosTrkChi2vNPoints_Pos;
       TH2F*  fMinosTrkChi2vNPoints_Neg;

        
        std::string fTrackMatchModuleLabel;
        std::string fMINOSModuleLabel;
        std::string fTracks_label;
        double fdBoundary; //distance from a boundary to be considered a track that "ends on a boundary"
        
    };
    
}



#endif // KINEMATICS_H
