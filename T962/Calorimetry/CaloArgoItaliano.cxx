//////////////////////////////////////////////////
//
// CaloArgoItaliano class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port echurch@fnal.gov
//  This algorithm is designed to perform the calorimetric reconstruction 
//  of the 3D reconstructed tracks
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "T962/Calorimetry/CaloArgoItaliano.h"
#include "Utilities/LArProperties.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"

#include "T962/T962_Objects/MINOS.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

//-------------------------------------------------
calo::CaloArgoItaliano::CaloArgoItaliano(fhicl::ParameterSet const& pset)
   : fTrackModuleLabel(pset.get< std::string >("TrackModuleLabel")      )
   , fnhits3D         (pset.get< int         >("SpacePointDimension")   )
   , fnhitsIND        (pset.get< int         >("HitInductionDimension") )
   , fnhitsCOL        (pset.get< int         >("HitCollectionDimension"))
   , fMINOSModuleLabel(pset.get< std::string >("MINOSModuleLabel")      )
   , fTrackMatchModuleLabel(pset.get< std::string >("TrackMatchModuleLabel"))

{
   //  produces< std::vector<calo::Calorimetry> >()
}

//-------------------------------------------------
calo::CaloArgoItaliano::~CaloArgoItaliano()
{
  
}

//-------------------------------------------------
void calo::CaloArgoItaliano::beginJob()
{
   // get access to the TFile service
   art::ServiceHandle<art::TFileService> tfs;

   /// ROOT tree for Calorimetry 
   ftree = tfs->make<TTree>("CaloTree","CaloTree");
   fdEdx_Coll = tfs->make<TH1F>("dEdx_Coll","dEdx Coll (MeV/cm)",1000,0.,10.); 
   fbirk = tfs->make<TH2F>("Birk_Coll_correction","Birk Coll MeV/cm-v-e/cm",100,0.0,5.0e5,20,0.0,50.0);
   fdEdx_Coll_vsXZangle = tfs->make<TH2F>("dEdx_Coll_vsXZ Angle","dEdx Coll vs. XZ Angle",1000,0.0,10.0,200,-1.0*TMath::Pi(),TMath::Pi());
   fdEdx_vs_Range = tfs->make<TH2F>("dEdx_vs_Range","dEdx vs. Range ",500,0.0,100.0,200,0.0,40.0);
   // dE/dx vs residual range plots for entering, contained, passing and escaping tracks (even if only for entering and contained tracks
   // the stopping point is know and the actual residual range can be calculated!) 
   fdEdx_vs_Range_ent  = tfs->make<TH2F>("dEdx_vs_Range_ent","dEdx vs. Range entering tracks",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_cont = tfs->make<TH2F>("dEdx_vs_Range_cont","dEdx vs. Range contained tracks",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_pass = tfs->make<TH2F>("dEdx_vs_Range_pass","dEdx vs. Range passing tracks",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_esc  = tfs->make<TH2F>("dEdx_vs_Range_esc","dEdx vs. Range escaping tracks",500,0.0,100.0,200,0.0,40.0);
   //    fdEdx_vs_Range_thr = tfs->make<TH2F>("dEdx_vs_Range_thr","dEdx vs. Range 4 > MeV/cm tracks",500,0.0,100.0,200,0.0,40.0);
   fKinetic_En = tfs->make<TH1F>("Kinetic_En","Kinetic Energy Deposited in LAr (MeV)",100,0.,500.); 

   fdEdx_vs_Range_MINOS_Pos = tfs->make<TH2F>("dEdx_vs_Range_MINOS_Pos","dEdx vs. Range : MINOS Pos",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_esc_MINOS_Pos  = tfs->make<TH2F>("dEdx_vs_Range_esc_MINOS_Pos","dEdx vs. Range escaping tracks: MINOS Pos",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_pass_MINOS_Pos  = tfs->make<TH2F>("dEdx_vs_Range_pass_MINOS_Pos","dEdx vs. Range passing tracks: MINOS Pos",500,0.0,100.0,200,0.0,40.0);
   fKinetic_En_MINOS_Pos = tfs->make<TH1F>("Kinetic_En_MINOS_Pos","Kinetic Energy Deposited in LAr (MeV): MINOS Pos",100,0.,500.);
   
   fdEdx_vs_Range_MINOS_Neg = tfs->make<TH2F>("dEdx_vs_Range_MINOS_Neg","dEdx vs. Range : MINOS Neg",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_esc_MINOS_Neg  = tfs->make<TH2F>("dEdx_vs_Range_esc_MINOS_Neg","dEdx vs. Range escaping tracks: MINOS Neg",500,0.0,100.0,200,0.0,40.0);
   fdEdx_vs_Range_pass_MINOS_Neg  = tfs->make<TH2F>("dEdx_vs_Range_pass_MINOS_Neg","dEdx vs. Range passing tracks: MINOS Neg",500,0.0,100.0,200,0.0,40.0);
   fKinetic_En_MINOS_Neg = tfs->make<TH1F>("Kinetic_En_MINOS_Neg","Kinetic Energy Deposited in LAr (MeV): MINOS Neg",100,0.,500.);

   fXHit = new double[fnhits3D];
   fYHit = new double[fnhits3D];
   fZHit = new double[fnhits3D];
   fwireIND = new int[fnhitsIND];
   ftimeIND = new double[fnhitsIND];
   fstimeIND = new double[fnhitsIND];
   fetimeIND = new double[fnhitsIND];
   fMIPsIND = new double[fnhitsIND];
   ftimeCOL = new double[fnhitsCOL];
   fstimeCOL = new double[fnhitsCOL];
   fetimeCOL = new double[fnhitsCOL];
   fMIPsCOL = new double[fnhitsCOL];
   fdEdxCOL = new double[fnhitsCOL];


   ftree->Branch("run", &frun, "run/I");
   ftree->Branch("event", &fevent, "event/I");
   ftree->Branch("itrack", &fntrack, "itrack/I");
   ftree->Branch("ntracks", &ftotTracks, "ntracks/I");
   ftree->Branch("TrkPitchI", &fTrkPitchI, "TrkPitchI/F");
   ftree->Branch("TrkPitchC", &fTrkPitchC, "TrkPitchC/F");
   ftree->Branch("XStart", &fXStart, "XStart/F");
   ftree->Branch("YStart", &fYStart, "YStart/F");
   ftree->Branch("ZStart", &fZStart, "ZStart/F");  
   ftree->Branch("XEnd", &fXEnd, "XEnd/F");  
   ftree->Branch("YEnd", &fYEnd, "YEnd/F");  
   ftree->Branch("ZEnd", &fZEnd, "ZEnd/F");  

   ftree->Branch("nhits3D",&fnhits3D,"nhits3D/I");
   ftree->Branch("XHit",fXHit,"XHit[nhits3D]/F");
   ftree->Branch("YHit",fYHit,"YHit[nhits3D]/F");
   ftree->Branch("ZHit",fZHit,"ZHit[nhits3D]/F");
   //  ftree->Branch("MIPs3D",&fMIPs3D,"MIPS3D[nhits3D]/F");

   ftree->Branch("nhitsIND",&fnhitsIND,"nhitsIND/I");
   ftree->Branch("wireIND",fwireIND,"wireIND[nhitsIND]/I");
   ftree->Branch("timeIND",ftimeIND,"timeIND[nhitsIND]/F");
   ftree->Branch("stimeIND",fstimeIND,"stimeIND[nhitsIND]/F");
   ftree->Branch("etimeIND",fetimeIND,"etimeIND[nhitsIND]/F");
   ftree->Branch("MIPsIND",fMIPsIND,"MIPsIND[nhitsIND]/F");

   ftree->Branch("nhitsCOL",&fnhitsCOL,"nhitsCOL/I");
   //  ftree->Branch("wireCOL",fwireCOL,"fwireCOL[nhitsCOL]/I");
   ftree->Branch("timeCOL",ftimeCOL,"timeCOL[nhitsCOL]/F");
   ftree->Branch("stimeCOL",fstimeCOL,"stimeCOL[nhitsCOL]/F");
   ftree->Branch("etimeCOL",fetimeCOL,"etimeCOL[nhitsCOL]/F");
   ftree->Branch("MIPsCOL",fMIPsCOL,"MIPsCOL[nhitsCOL]/F");
   ftree->Branch("dEdxCOL",fdEdxCOL,"dEdxCOL[nhitsCOL]/F");

}

//------------------------------------------------------------------------------------//
void calo::CaloArgoItaliano::analyze(const art::Event& evt)
{ 
  
   art::Handle< std::vector<recob::Track> > trackListHandle;
   evt.getByLabel(fTrackModuleLabel,trackListHandle);

   art::PtrVector<recob::Track> tracklist;  
    
   for (unsigned int ii = 0; ii <  trackListHandle->size(); ++ii)
   {
      art::Ptr<recob::Track> track(trackListHandle,ii);
      tracklist.push_back(track);
   }

   // Conversion factors
   const double calFactor = 7.54;   //   in ADC/fC
   const double eCharge   = 1.6e-4; //   1e = 1.6 *10^-4 fC


   // Get Geometry
   art::ServiceHandle<geo::Geometry> geom;

   //TPC dimensions (Local reference frame) To be checked
   TPCsize[0] = (geom->DetHalfWidth()-5.)*2.;
   TPCsize[1] = (geom->DetHalfHeight()-0.5715)*2.;
   TPCsize[2] = (geom->DetLength()-5.)*1.;

   frun = evt.id().run();
   fevent = evt.id().event();

   fntrack = 0; 
   fnhits3D = 0; fnhitsIND = 0; fnhitsCOL = 0;

   art::Handle< std::vector<t962::MINOS> > MinosTrackHandle;
   evt.getByLabel(fMINOSModuleLabel,MinosTrackHandle);
   
   art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
   evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);
  
   for(art::PtrVector<recob::Track>::const_iterator trkIter = tracklist.begin(); trkIter != tracklist.end();  trkIter++){   

      
      bool match_positive =false;
      bool match_negative = false;
      for(unsigned int j=0; j<MinosTrackHandle->size();++j)
      {
         art::Ptr<t962::MINOS> minostrack(MinosTrackHandle,j);
         for (unsigned int k = 0; k < trackmatchListHandle->size(); k++){
            art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,k);
            if( (trackmatchHolder->fMINOStrackid == minostrack->ftrkIndex) &&
                (trackmatchHolder->fArgoNeuTtrackid == (*trkIter)->ID()) &&
                minostrack->ftrkmom>0.0){
               if(minostrack->fcharge==1) match_positive = true;
               if(minostrack->fcharge==-1) match_negative = true;
            }
         }
      }

      
      
     
             
      fnhits3D = ((*trkIter)->SpacePoints()).size();

      //store track directional cosines
      double trackCosStart[3]={0.,0.,0.};
      double trackCosEnd[3]={0.,0.,0.};
      (*trkIter)->Direction(trackCosStart,trackCosEnd);
        
      bool startsonboundary = BeginsOnBoundary(*trkIter);
      bool endsonboundary   = EndsOnBoundary(*trkIter);

      int containment = 0;
      if(!startsonboundary && !endsonboundary ) containment = 0;  // contained track
      if(!startsonboundary && endsonboundary )  containment = 1;  // escaping track
      if(startsonboundary && !endsonboundary )  containment = 2;   // entering track
      if(startsonboundary && endsonboundary )   containment = 3;  // passing track 
 
      std::cout<<"|-* ArgoNeuT track #"<<fntrack;
      switch (containment) {
         case 0:
            std::cout<<", contained"<<std::endl;
            break;
         case 1:
            std::cout<<", escaping"<<std::endl;
            break;
         case 2:
            std::cout<<", entering"<<std::endl;
            break;
         case 3:
            std::cout<<", passing"<<std::endl;
            break;
         default:
            std::cout<<", ??"<<std::endl;
      }
  
      //      fMIPs3D = new double[fnhits3D];

      //recover the Induction (geo::kU) cluster list
      art::PtrVector<recob::Cluster> kUclusterlist = (*trkIter)->Clusters(geo::kU);
    
      int npI = 0;
    
      // Some variables for the hit
      unsigned int channel;//channel number
      float time;          //hit time at maximum
      float stime;         //hit start time 
      float etime;         //hit end time 
      unsigned int tpc;    //hit tpc number 
      unsigned int wire;   //hit wire number 
      unsigned int plane;  //hit plane number
    
      // loop over all (geo::kU) clusters of the track
      for(art::PtrVector<recob::Cluster>::const_iterator clIter = kUclusterlist.begin(); 
          clIter != kUclusterlist.end();  
          clIter++) {   
      
         if((*clIter)->View() == geo::kU) fTrkPitchI = (*trkIter)->ProjectedLength(geo::kU);
                  
         art::PtrVector<recob::Hit> hitlist;
         hitlist = (*clIter)->Hits(geo::kU);
      
         fnhitsIND = npI+hitlist.size();
      
         //loop over cluster hits
         for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
            //recover the Hit
            time = (*hitIter)->PeakTime() ;
            stime = (*hitIter)->StartTime() ;
            etime = (*hitIter)->EndTime();            
            art::Ptr <recob::Wire> theWire = (*hitIter)->Wire();
            channel = theWire->RawDigit()->Channel();
            geom->ChannelToWire(channel,tpc,plane,wire);
	
            double MIPs = (*hitIter)->Charge(true);
	
            fwireIND[npI] = (int)wire;
            ftimeIND[npI] = (double)time;
            fstimeIND[npI] = (double)stime;
            fetimeIND[npI] = (double)etime;
            fMIPsIND[npI] = (double)MIPs;
	
            npI++;
         }// end of loop over the (geo::kU) cluster hits
	  
      }// end of loop over all I clusters of the track
      
      int npC = 0;
    
      //recover the Collection (geo::kV) cluster list
      art::PtrVector<recob::Cluster> kVclusterlist = (*trkIter)->Clusters(geo::kV);
    
      // loop over all (geo::KV) clusters of the track
      for(art::PtrVector<recob::Cluster>::const_iterator clIter = kVclusterlist.begin(); 
          clIter != kVclusterlist.end();  
          clIter++) { 

         if((*clIter)->View() == geo::kV) fTrkPitchC = (*trkIter)->ProjectedLength(geo::kV);  

         // Some variables for the hit
         float time;         //hit time at maximum
         float stime;        //hit start time 
         float etime;        //hit end time 
         art::PtrVector<recob::Hit> hitlist;
         hitlist = (*clIter)->Hits(geo::kV);
      
         fnhitsCOL = npC + hitlist.size();
      
         double dEdx_Coll_Tot = 0.;
         double Kin_En = 0.; 
         //loop over cluster hits
         for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
            //recover the Hit
            time = (*hitIter)->PeakTime() ;
            stime = (*hitIter)->StartTime() ;
            etime = (*hitIter)->EndTime();            
            art::Ptr<recob::Wire> theWire = (*hitIter)->Wire();
            channel = theWire->RawDigit()->Channel();
            geom->ChannelToWire(channel,tpc,plane,wire);
	      
            double MIPs   = (*hitIter)->Charge(true);            // in ADC
            double dQdx   = MIPs/fTrkPitchC;           // in ADC/cm
            double dQdx_e = dQdx/(calFactor*eCharge);  // in e/cm
            if(evt.isRealData()) dQdx_e *= LifetimeCorrection(frun,time);   // in e/cm
            double dEdx    = BirksCorrection(dQdx_e);         // in MeV/cm
            //       std::cout<<"MIPs (ADC) "<<MIPs<<" dQdx (ADC/cm) "<<dQdx<<" corr "<<LifetimeCorrection(frun,time)<<" dQdx_e (e/cm)"<<dQdx_e<<" dEdx (MeV/cm) "<<dEdx<<std::endl;
            fdEdx_Coll->Fill(dEdx);
            fbirk->Fill(dQdx_e,dEdx);
            fdEdx_Coll_vsXZangle->Fill(dEdx,TMath::ATan(trackCosStart[0]/trackCosStart[2]));
	      
            dEdx_Coll_Tot = dEdx_Coll_Tot + dEdx;
            Kin_En = Kin_En + dEdx * fTrkPitchC;
	      
            ftimeCOL[npC] = (double)time;
            fstimeCOL[npC] = (double)stime;
            fetimeCOL[npC] = (double)etime;
            fMIPsCOL[npC] =  (double)MIPs;
            fdEdxCOL[npC] = (double)dEdx;
            npC++;
	
            // dE/dx vs Range plots
            for(std::vector<recob::SpacePoint>::const_iterator spIter = ((*trkIter)->SpacePoints()).begin(); 
                spIter < ((*trkIter)->SpacePoints()).end();  
                ++spIter) { 
               art::PtrVector<recob::Hit> sphits = (*spIter).Hits(geo::kV);
               if(sphits.size()==0) continue;
               if((*hitIter)->Channel()!=sphits[0]->Channel()) continue;
               if((*hitIter)->PeakTime()!=sphits[0]->PeakTime()) continue;//protect against multiple hits on the same wire in the cluster.
               const double *xyz = new double[3];
               xyz = (*spIter).XYZ();
               std::vector<double> larStart, larEnd;
               (*trkIter)->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)
               double range = sqrt((larEnd[0]-xyz[0])*(larEnd[0]-xyz[0]) +
                                   (larEnd[1]-xyz[1])*(larEnd[1]-xyz[1]) +
                                   (larEnd[2]-xyz[2])*(larEnd[2]-xyz[2]));

               fdEdx_vs_Range->Fill(range,dEdx);
               if(containment == 0) fdEdx_vs_Range_cont->Fill(range,dEdx); 
               if(containment == 1) fdEdx_vs_Range_esc->Fill(range,dEdx);   
               if(containment == 2) fdEdx_vs_Range_ent->Fill(range,dEdx);   
               if(containment == 3) fdEdx_vs_Range_pass->Fill(range,dEdx); 

               if(match_positive){
                  fdEdx_vs_Range_MINOS_Pos->Fill(range,dEdx);
                  if(containment == 1) fdEdx_vs_Range_esc_MINOS_Pos->Fill(range,dEdx);    
                  if(containment == 3) fdEdx_vs_Range_pass_MINOS_Pos->Fill(range,dEdx); 
               }
               if(match_negative){
                  fdEdx_vs_Range_MINOS_Neg->Fill(range,dEdx);
                  if(containment == 1) fdEdx_vs_Range_esc_MINOS_Neg->Fill(range,dEdx);    
                  if(containment == 3) fdEdx_vs_Range_pass_MINOS_Neg->Fill(range,dEdx); 
               }
                     
               //npC++;
               
            }
         }// end of loop over the (geo::kV) cluster hits
      
         double dEdx_Coll_Mean = dEdx_Coll_Tot/npC;
            
         fKinetic_En->Fill(Kin_En);
         if(match_positive) fKinetic_En_MINOS_Pos->Fill(Kin_En);
         if(match_negative) fKinetic_En_MINOS_Neg->Fill(Kin_En);


//             if(dEdx_Coll_Mean > 4. ) {
// 	      //loop over cluster hits
// 	      for(art::PtrVectorItr<recob::Hit> hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
//                 //recover the Hit
//                 time = (*hitIter)->PeakTime() ;
//                 stime = (*hitIter)->StartTime() ;
//                 etime = (*hitIter)->EndTime();            
//                 art::Ptr<recob::Wire> theWire = (*hitIter)->Wire();
//                 channel = theWire->RawDigit()->Channel();
//                 geom->ChannelToWire(channel,tpc,plane,wire);
	
//                 double MIPs   = (*hitIter)->Charge(true);            // in ADC
//                 double dQdx   = MIPs/fTrkPitchC;           // in ADC/cm
//                 double dQdx_e = dQdx/(calFactor*eCharge);  // in e/cm
//                 if(evt.isRealData()) dQdx_e *= LifetimeCorrection(frun,time);   // in e/cm
//                 double dEdx    = BirksCorrection(dQdx_e);         // in MeV/cm
//                 //	      std::cout<<"MIPs (ADC) "<<MIPs<<" dQdx (ADC/cm) "<<dQdx<<" corr "<<LifetimeCorrection(frun,time)<<" dQdx_e (e/cm)"<<dQdx_e<<" dEdx (MeV/cm) "<<dEdx<<std::endl;

// 		if(!endsonboundary) {
// 		  for(std::vector<recob::SpacePoint>::const_iterator spIter = ((*trkIter)->SpacePoints()).begin(); 
// 		      spIter < ((*trkIter)->SpacePoints()).end();  
// 		      ++spIter) { 
// 		    art::PtrVector<recob::Hit> sphits = (*spIter).Hits(geo::kV);
// 		    if(sphits.size()==0) continue;
// 		    if((*hitIter)->Channel()!=sphits[0]->Channel()) continue;
// 		    const double *xyz = new double[3];
// 		    xyz = (*spIter).XYZ();
// 		    std::vector<double> larStart, larEnd;
// 		    (*trkIter)->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)
// 		    double range = sqrt((larEnd[0]-xyz[0])*(larEnd[0]-xyz[0]) +
// 					(larEnd[1]-xyz[1])*(larEnd[1]-xyz[1]) +
// 					(larEnd[2]-xyz[2])*(larEnd[2]-xyz[2]));
// 		    fdEdx_vs_Range_thr->Fill(range,dEdx);
// 		  }
// 		}
// 	      } // end of loop over the (geo::kV) cluster hits
// 	    } 

         std::cout<<"  |-* Collection View Calorimetric Reco"<<std::endl;
         std::cout<<"    |-* Hits="<<fnhitsCOL<<std::endl;
         std::cout<<"    |-* <dE/dx>="<<dEdx_Coll_Mean<<" MeV/cm"<<std::endl;
         std::cout<<"    |-* Kinetic Energy deposited in LAr="<<Kin_En<<" MeV"<<std::endl;

      }// end of loop over all C clusters of the track
            
    
      ftotTracks = (int ) (tracklist.size());
      ftree->Fill();

      fntrack++;
    
   }// end of loop over all tracks created by  Track3Dreco

   return;
}

//------------------------------------------------------------------------------------//
double calo::CaloArgoItaliano::LifetimeCorrection(int Run, float time){
   int run_num = Run;
   float t = time;
   art::ServiceHandle<util::LArProperties> LArProp;

   double timetick = 0.198;    //time sample in us
   double presamplings = 60.;
   double plane_pitch = 0.4;   //wire plane pitch in cm 
   double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
   double Efield_IC = 0.89;     // Electric Field between Induction and Collection planes in kV/cm
   double Temperature = 87.6;  // LAr Temperature in K
   double driftvelocity = 0.157;                                      // from the measurement on passing muons
   // DriftVelocity() in Utilities. EC, 6-Dec-2010.
   double driftvelocity_SI = LArProp->DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
   double driftvelocity_IC = LArProp->DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
   double timepitch = driftvelocity*timetick;                         //time sample (cm) 
   double tSI = plane_pitch/driftvelocity_SI/timetick;                //drift time between Shield and Collection planes (time samples)
   double tIC = plane_pitch/driftvelocity_IC/timetick;                //drift time between Induction and Collection planes (time samples)

   t -= presamplings;
   t -= (tSI+tIC);   // Collection
   double t_cm = (Double_t)(t * timepitch);
   double tau_cm = t_cm;
   double correction;

   // needs a Data Base for the extraction of lifetime..
   if(run_num==609) tau_cm = 732.*driftvelocity; //  (in cm)
   if(run_num==618) tau_cm = 732.*driftvelocity; //  (in cm)
   if(run_num==620) tau_cm = 710.*driftvelocity; //  (in cm)  
   if(run_num==621) tau_cm = 734.*driftvelocity; //  (in cm)  
   if(run_num==627) tau_cm = 757.*driftvelocity; //  (in cm)  
   if(run_num==628) tau_cm = 765.*driftvelocity; //  (in cm)  
   if(run_num==629) tau_cm = 789.*driftvelocity; //  (in cm)  
   if(run_num==632) tau_cm = 789.*driftvelocity; //  (in cm)  
   if(run_num==633) tau_cm = 737.*driftvelocity; //  (in cm)  
   if(run_num==634) tau_cm = 744.*driftvelocity; //  (in cm)  
   if(run_num==635) tau_cm = 764.*driftvelocity; //  (in cm) 
   if(run_num==640) tau_cm = 769.*driftvelocity; //  (in cm)  
   if(run_num==644) tau_cm = 763.*driftvelocity; //  (in cm)  
   if(run_num==648) tau_cm = 772.*driftvelocity; //  (in cm)   
   if(run_num==649) tau_cm = 764.*driftvelocity; //  (in cm)  
   if(run_num==650) tau_cm = 735.*driftvelocity; //  (in cm) 
     
   correction = exp(t_cm/tau_cm);
   return correction;
}
//------------------------------------------------------------------------------------//
double calo::CaloArgoItaliano::BirksCorrection(double dQdx_e){

   /// Correction for charge quenching using parameterization from 
   /// S.Amoruso et al., NIM A 523 (2004) 275

   double dQdx = dQdx_e;
   double dEdx;
    
   float A3t = 0.800; 
   float K3t = 0.0486; // in KV/cm*(g/cm^2)/MeV
   float rho = 1.388; // LAr density in g/cm^3
   double Wion = 23.6e-6 ; //   23.6 eV = 1e, Wion in MeV/e
   double Efield = 0.485;  // Electric Field in the drift region in KV/cm
   K3t = K3t/rho; // KV/MeV
   dEdx = dQdx/(A3t/Wion-K3t/Efield*dQdx); //MeV/cm 
  
   return dEdx;
}

//------------------------------------------------------------------------------------//
void calo::CaloArgoItaliano::ReadCaloTree(){
   // This method is an example of how to read 
   // the dynamic vectors in the ROOT Tree

   ftree->Scan("run:event:TrkPitchI:TrkPitchC:XStart:nhits3D:nhitsIND:nhitsCOL");
   int nentries=(int)ftree->GetEntries();
   std::cout<<"nentries  "<<nentries<<std::endl;

   std::vector<double> *MIPsCOL = 0;
   TBranch *bMIPsCOL = 0;
   ftree->SetBranchAddress("MIPsCOL",&MIPsCOL,&bMIPsCOL);

   std::vector<double> *MIPsIND = 0;
   TBranch *bMIPsIND = 0;
   ftree->SetBranchAddress("MIPsIND",&MIPsIND,&bMIPsIND);
   for (int i = 0; i < nentries; i++) {
      std::cout<<" entry "<<i<<std::endl;
      long int tentry = ftree->LoadTree(i);

      bMIPsCOL->GetEntry(tentry);
      std::cout<<"# of hits COLL "<<MIPsCOL->size()<<std::endl;
      for (unsigned int j = 0; j < MIPsCOL->size(); ++j) {
         std::cout<<" Coll "<<MIPsCOL->at(j)<<std::endl;          
      }
      bMIPsIND->GetEntry(tentry);
      std::cout<<"# of hits IND "<<MIPsIND->size()<<std::endl;
      for (unsigned int j = 0; j < MIPsIND->size(); ++j) {
         std::cout<<" Ind "<<MIPsIND->at(j)<<std::endl;          
      }
   }
}

//--------------------------------------------------
bool calo::CaloArgoItaliano::BeginsOnBoundary(art::Ptr<recob::Track> lar_track)
{
   double fdBoundary = 1.5;
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
bool calo::CaloArgoItaliano::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
{
   double fdBoundary = 1.5;
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

//-------------------------------
