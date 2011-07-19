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

        TH1D* fMinosErange_Pos;
        TH1D* fMinosErange_Neg;
        TH1D* fMinosMom_Pos;
        TH1D* fMinosMom_Neg;

        
        std::string fTrackMatchModuleLabel;
        std::string fMINOSModuleLabel;
        std::string fTracks_label;
        double fdBoundary; //distance from a boundary to be considered a track that "ends on a boundary"
        
    };
    
}



#endif // KINEMATICS_H
