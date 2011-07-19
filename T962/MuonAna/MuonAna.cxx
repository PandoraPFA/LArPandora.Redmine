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
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "cetlib/exception.h"

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
#include "T962/T962_Objects/MINOSTrackMatch.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

namespace muons {
//-----------------------------------------------------------------------------
    MuonAna::MuonAna(fhicl::ParameterSet const& pset) :
        fMINOSModuleLabel      (pset.get< std::string > ("MINOSModuleLabel")),
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
        fVertAngle = tfs->make<TH1F>("fVertAngle","VertAngle", 200,-1.0*TMath::Pi(),TMath::Pi());
        fHorizAngle = tfs->make<TH1F>("fHorizAngle","HorizAngle", 200,-1.0*TMath::Pi(),TMath::Pi());

        fStartXvsStartY = tfs->make<TH2F>("fStartXvsStartY","Start X vs. Y", 220,-5.0,50.0,200,-25.0,25.0);
        fStartZvsStartX = tfs->make<TH2F>("fStartZvsStartX","Start Z vs. X", 400,-5.0,95.0,220,-5.0,50.0);
        fStartZvsStartY = tfs->make<TH2F>("fStartZvsStartY","Start Z vs. Y", 400,-5.0,95.0,200,-25.0,25.0);
        fEndXvsEndY = tfs->make<TH2F>("fEndXvsEndY","End X vs. Y", 220,-5.0,50.0,200,-25.0,25.0);
        fEndZvsEndX = tfs->make<TH2F>("fEndZvsEndX","End Z vs. X", 400,-5.0,95.0,220,-5.0,50.0);
        fEndZvsEndY = tfs->make<TH2F>("fEndZvsEndY","End Z vs. Y", 400,-5.0,95.0,200,-25.0,25.0);

        fMinosErange_Pos = tfs->make<TH1D>("fMinosErange_Pos","MINOS + Charge Tracks: Erange",1000,0.0,10000.0);
        fMinosErange_Neg = tfs->make<TH1D>("fMinosErange_Neg","MINOS - Charge Tracks: Erange",1000,0.0,10000.0);
        fMinosMom_Pos = tfs->make<TH1D>("fMinosMom_Pos","MINOS + Charge Tracks: Momentum",1000,0.0,10000.0);
        fMinosMom_Neg = tfs->make<TH1D>("fMinosMom_Neg","MINOS - Charge Tracks: Momentum",1000,0.0,10000.0);

            
        
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

        
        if(LarTrackHandle->size()>0){
            for(unsigned int i=0; i<LarTrackHandle->size();++i){
                
                art::Ptr<recob::Track> lartrack(LarTrackHandle,i);
              
                
                bool startsonboundary = BeginsOnBoundary(lartrack);
                bool endsonboundary   = EndsOnBoundary(lartrack);

                if(startsonboundary && endsonboundary){//throughgoing track...make plots
                   //std::cout << " MuonAna : " << *lartrack << std::endl;
                    lartrack->Extent(larStart,larEnd);
                    lartrack->Direction(trackCosStart,trackCosEnd);

                    double tracklength = sqrt((larStart[0]-larEnd[0])*(larStart[0]-larEnd[0]) +
                                              (larStart[1]-larEnd[1])*(larStart[1]-larEnd[1]) +
                                              (larStart[2]-larEnd[2])*(larStart[2]-larEnd[2]));

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
                    fVertAngle->Fill(TMath::ATan(trackCosStart[1]/trackCosStart[2]));
                    fHorizAngle->Fill(TMath::ATan(trackCosStart[0]/trackCosStart[2]));
                    fStartXvsStartY->Fill(larStart[0],larStart[1]);
                    fStartZvsStartX->Fill(larStart[2],larStart[0]);
                    fStartZvsStartY->Fill(larStart[2],larStart[1]);
                    fEndXvsEndY->Fill(larEnd[0],larEnd[1]);
                    fEndZvsEndX->Fill(larEnd[2],larEnd[0]);
                    fEndZvsEndY->Fill(larEnd[2],larEnd[1]);

                    art::Handle< std::vector<t962::MINOS> > MinosTrackHandle;
                    evt.getByLabel(fMINOSModuleLabel,MinosTrackHandle);

                    art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
                    evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);


                    for(unsigned int j=0; j<MinosTrackHandle->size();++j)
                    {
                        art::Ptr<t962::MINOS> minostrack(MinosTrackHandle,j);

                        bool match = false;
                        for (unsigned int k = 0; k < trackmatchListHandle->size(); k++){
                            art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,k);
                            if( (trackmatchHolder->fMINOStrackid == minostrack->ftrkIndex) &&
                                (trackmatchHolder->fArgoNeuTtrackid == lartrack->ID())){
                                if(minostrack->fcharge==1){
                                    fMinosErange_Pos->Fill(1000.0*minostrack->ftrkErange);
                                    fMinosMom_Pos->Fill(1000.0*minostrack->ftrkmom);
                                }
                                if(minostrack->fcharge==-1){
                                    fMinosErange_Neg->Fill(1000.0*minostrack->ftrkErange);
                                    fMinosMom_Neg->Fill(1000.0*minostrack->ftrkmom);
                                }
                            }//found matching tracks
                        }//trackmatch loop
                    }//MINOS loop
                  
                    
                    
                    
                }
                
            }
        }

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
