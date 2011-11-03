#ifndef CALO_H
#define CALO_H

#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h" // include the proper bit of the framework

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include "RecoBase/Track.h"

///calorimetry
namespace calo {
   
    class CaloArgoItaliano : public art::EDAnalyzer {
    
    public:
    
        explicit CaloArgoItaliano(fhicl::ParameterSet const& pset); 
        virtual ~CaloArgoItaliano();
    
        void beginJob(); 
        //    void endJob();

        void analyze(const art::Event& evt);

    private:
        
        double LifetimeCorrection(int Run, float time);
        double BirksCorrection(double dQdx_e);
	void   ReadCaloTree();

        bool BeginsOnBoundary(art::Ptr<recob::Track> lar_track);
        bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);
    
        std::string     fTrackModuleLabel;
        std::string fTrackMatchModuleLabel;
        std::string fMINOSModuleLabel;

        TTree *ftree;
        TH1F *fdEdx_Coll;
        TH2F *fbirk;
        TH2F *fdEdx_Coll_vsXZangle;
	TH2F *fdEdx_vs_Range;
	TH2F *fdEdx_vs_Range_ent;
	TH2F *fdEdx_vs_Range_cont;
	TH2F *fdEdx_vs_Range_pass;
	TH2F *fdEdx_vs_Range_esc;
	TH2F *fdEdx_vs_Range_thr;
       TH1F *fKinetic_En;

       TH2F *fdEdx_vs_Range_MINOS_Pos;
       TH2F *fdEdx_vs_Range_esc_MINOS_Pos;
       TH2F *fdEdx_vs_Range_pass_MINOS_Pos;
       TH1F* fKinetic_En_MINOS_Pos;
       TH1F* fKinetic_En_esc_MINOS_Pos;
       TH1F* fKinetic_En_pass_MINOS_Pos;
       TH2F *fdEdx_vs_Range_MINOS_Neg;
       TH2F *fdEdx_vs_Range_esc_MINOS_Neg;
       TH2F *fdEdx_vs_Range_pass_MINOS_Neg;
       TH1F* fKinetic_En_MINOS_Neg;
       TH1F* fKinetic_En_esc_MINOS_Neg;
       TH1F* fKinetic_En_pass_MINOS_Neg;
 
        int frun;           //Run 
        int fevent;         //Event
        int fntrack;        //track number
        int ftotTracks;        //track number
        double fTrkPitchI;  //
        double fTrkPitchC;
        double fXStart;
        double fYStart;
        double fZStart;
        double fXEnd;
        double fYEnd;
        double fZEnd;

        int fnhits3D; 
        double* fXHit;
        double* fYHit;
        double* fZHit;
        double* fMIPs3D;   
 
        int fnhitsIND; 
        int* fwireIND;
        double* ftimeIND;
        double* fstimeIND;
        double* fetimeIND;
        double* fMIPsIND;

        int fnhitsCOL; 
        int* fwireCOL;
        double* ftimeCOL;
        double* fstimeCOL;
        double* fetimeCOL;
        double* fMIPsCOL;
        double* fdEdxCOL;

 
        double TPCsize[3];

    protected: 
    
  
    }; // class CaloArgoItaliano

}

#endif // CALO_H
