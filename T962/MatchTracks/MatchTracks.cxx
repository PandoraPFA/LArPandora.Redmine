////////////////////////////////////////////////////////////////////////
//
//   MatchTracks Class
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/AssociationUtil.h"

#include "TH1D.h"
#include "TH2D.h"


#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <dirent.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>


#include "T962/MatchTracks/MatchTracks.h"

namespace match{

   //-------------------------------------------------
   MatchTracks::MatchTracks(fhicl::ParameterSet const& pset) 
   {
      this->reconfigure(pset);

      produces< art::Assns<recob::Track, t962::MINOS> >();
   }

   //-------------------------------------------------
   MatchTracks::~MatchTracks()
   {
   
   }

   //-------------------------------------------------
   void MatchTracks::reconfigure(fhicl::ParameterSet const& pset)
   {
      fLarTracks_label = pset.get< std::string >("lartracks");
      fMinosTracks_label = pset.get< std::string >("minostracks");
      fdZ = pset.get< double>("dZ");
      fdXY = pset.get< double>("dXY");
      fdCosx = pset.get< double >("dCosx");
      fdCosy = pset.get< double >("dCosy");
      fdCosz = pset.get< double >("dCosz");
      fdBoundary = pset.get< double >("dBoundary");
   }

   //-------------------------------------------------
   void MatchTracks::beginJob()
   {
      //get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;

      fDiffDirCosX = tfs->make<TH1D>("fDiffDirCosX","Diff: DirCosX", 100,-2.0,2.0);
      fDiffDirCosY = tfs->make<TH1D>("fDiffDirCosY","Diff: DirCosY", 100,-2.0,2.0);
      fDiffDirCosZ = tfs->make<TH1D>("fDiffDirCosZ","Diff: DirCosZ", 100,-2.0,2.0);
      fDiffR       = tfs->make<TH1D>("fDiffR","Diff: Radial XY",200,0.0,200.0);
      fDiffTotal       = tfs->make<TH1D>("fDiffTotal","Diff: Radial XY/cos(theta)",200,0.0,200.0);
      fDiffXvDiffY = tfs->make<TH2D>("fDiffXvDiffY","Diff: X v Y",200,-100.0,100.0,200,-100.0,100.0);

      fT962_DirCosX   = tfs->make<TH1D>("fT962_DirCosX","T962 Track DirCosX",100, -1.0,1.0); 
      fT962_DirCosY   = tfs->make<TH1D>("fT962_DirCosY","T962 Track DirCosY",100, -1.0,1.0); 
      fT962_DirCosZ   = tfs->make<TH1D>("fT962_DirCosZ","T962 Track DirCosZ",100, -1.0,1.0); 
      fT962_StartZY =   tfs->make<TH2D>("fT962_StartZY","T962 Track Start ZY",100,-5.0,95.0,50,-25.0,25.0);
      fT962_EndZY =     tfs->make<TH2D>("fT962_EndZY","T962 Track End ZY",100,-5.0,95.0,50,-25.0,25.0);
      fT962_StartXY =   tfs->make<TH2D>("fT962_StartXY","T962 Track Start XY",60,-5.0,55.0,50,-25.0,25.0);
      fT962_EndXY =     tfs->make<TH2D>("fT962_EndXY","T962 Track End XY",60,-5.0,55.0,50,-25.0,25.0);
      fT962_StartZX =   tfs->make<TH2D>("fT962_StartZX","T962 Track Start ZX",100,-5.0,95.0,60,-5.0,55.0);
      fT962_EndZX =     tfs->make<TH2D>("fT962_EndZX","T962 Track End ZX",100,-5.0,95.0,60,-5.0,55.0);
      fT962_Ntracks =   tfs->make<TH1D>("fT962_Ntracks","T962 # Tracks",20,0.0,20.0);
      fMINOS_Ntracks =   tfs->make<TH1D>("fMINOS_Ntracks","MINOS # Tracks",20,0.0,20.0);
      fMinos_DirCosX   = tfs->make<TH1D>("fMinos_DirCosX","MINOS Track DirCosX",100, -1.0,1.0); 
      fMinos_DirCosY   = tfs->make<TH1D>("fMinos_DirCosY","MINOS Track DirCosY",100, -1.0,1.0); 
      fMinos_DirCosZ   = tfs->make<TH1D>("fMinos_DirCosZ","MINOS Track DirCosZ",100, -1.0,1.0); 
      fMinos_XY     = tfs->make<TH2D>("fMinos_XY","MINOS Start XY",600,-300.0,300.0,600,-300.0,300.0);
      fMinosErange_Pos = tfs->make<TH1D>("fMinosErange_Pos","MINOS + Charge Tracks: Erange",1000,0.0,10000.0);
      fMinosErange_Neg = tfs->make<TH1D>("fMinosErange_Neg","MINOS - Charge Tracks: Erange",1000,0.0,10000.0);
      fMinosMom_Pos = tfs->make<TH1D>("fMinosMom_Pos","MINOS + Charge Tracks: Momentum",1000,0.0,10000.0);
      fMinosMom_Neg = tfs->make<TH1D>("fMinosMom_Neg","MINOS - Charge Tracks: Momentum",1000,0.0,10000.0);

      fDiffXvD = tfs->make<TH2D>("fDiffXvD","Diff: X Pred. - MINOS x",100,-200.0,200.0,100,90.0,190.0);
      fDiffYvD = tfs->make<TH2D>("fDiffYvD","Diff: Y Pred. - MINOS y",100,-200.0,200.0,100,90.0,190.0);
   }

   //-------------------------------------------------
   void MatchTracks::produce(art::Event& evt)
   {
      std::unique_ptr< art::Assns<recob::Track, t962::MINOS> > assn(new art::Assns<recob::Track, t962::MINOS>);
   
      art::Handle< std::vector<recob::Track> > LarTrackHandle;
      if (!evt.getByLabel(fLarTracks_label,LarTrackHandle)) return;

      art::Handle< std::vector<t962::MINOS> > MinosTrackHandle;
      if (!evt.getByLabel(fMinosTracks_label,MinosTrackHandle)) return;

      mf::LogVerbatim("MatchTracks") << std::setfill('*') << std::setw(175) << "*" << std::setfill(' ') << "\n"
                                 << " #T962 Tracks = " << LarTrackHandle->size() 
                                 << " #MINOS Tracks = " << MinosTrackHandle->size() ;

      fT962_Ntracks->Fill(LarTrackHandle->size());
      fMINOS_Ntracks->Fill(MinosTrackHandle->size());
      int exiting = 0;
      for(unsigned int i=0; i<LarTrackHandle->size();++i){
         art::Ptr<recob::Track> lartrack(LarTrackHandle,i);
         if(EndsOnBoundary(lartrack)) ++exiting;
         mf::LogVerbatim("MatchTracks") << "T962 " << *lartrack ;
      }
      mf::LogVerbatim("MatchTracks") << exiting << " Exiting T962 tracks." ;

      for(unsigned int j=0; j<MinosTrackHandle->size();++j){
         art::Ptr<t962::MINOS> minostrack(MinosTrackHandle,j);
         mf::LogVerbatim("MatchTracks") << *minostrack ;
      }
      
      //the following makes sure that only 1 argoneut track is assigned to one minos track and that the match is the strongest (in terms of projected radial difference and angle between the tracks) among the candidate matches       
      for(unsigned int i=0; i<LarTrackHandle->size();++i)
      {
         double totaldiff2=1000000.;
         double rdiff_best = -999.0;
         double xdiff_best = -999.0;
         double ydiff_best = -999.0;
         unsigned int matchnumber = 999;
         art::Ptr<recob::Track> lartrack(LarTrackHandle,i);
         if(!EndsOnBoundary(lartrack)) continue;//track doesn't leave TPC

         std::vector<double> larStart, larEnd;
         lartrack->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)            
         double lardirectionStart[3];
         double lardirectionEnd[3];
         lartrack->Direction(lardirectionStart,lardirectionEnd); 
         
         for(unsigned int j=0; j<MinosTrackHandle->size();++j)
         {
            art::Ptr<t962::MINOS> minostrack(MinosTrackHandle,j);            
            if((100.0*minostrack->ftrkVtxZ)>fdZ) continue;//MINOS track too far into detector

            //some plots to check alignment
            for(int n = 0; n<100;++n){
               //double D=(90*0.5)+(42.4*2.54)-5.588;
               double D = 90.0 + n*1.0;
               double dz = D - larEnd[2]+(100.0 * minostrack->ftrkVtxZ);//z-difference between end of T962 track and
               
               double l = dz/(lardirectionEnd[2]);//3-d distance between end of T962 track and begin of MINOS track
               
               double x_pred = l*lardirectionEnd[0]+larEnd[0];//predicted x-pos. of T962 track at z-position equal to
               //start of MINOS track
               double y_pred = l*lardirectionEnd[1]+larEnd[1];//predicted y-pos. of T962 track at z-position equal to
               //start of MINOS track
               
               fDiffXvD->Fill((100.0*minostrack->ftrkVtxX - x_pred),D);
               fDiffYvD->Fill((100.0*minostrack->ftrkVtxY - y_pred),D);
            }

                  
            double xdiff,ydiff,rdiff,totaldiff;   
            bool match = Compare(lartrack,minostrack,xdiff,ydiff,rdiff,totaldiff,totaldiff2);  
            
            if(match){
               matchnumber = j;//keep track of best match found so far
               rdiff_best = rdiff;
               xdiff_best = xdiff;
               ydiff_best = ydiff;
            }

         }//loop over MINOS tracks

         
         //Check if best match was found, then associate it to Track
         if(matchnumber!=999){
            art::Ptr<t962::MINOS> minostrack(MinosTrackHandle,matchnumber);
            mf::LogVerbatim("MatchTracks") << " Match! T962 Track #" << lartrack->ID() 
                                       << " and MINOS Track #" << minostrack->ftrkIndex;

            fDiffR->Fill(rdiff_best);
            fDiffTotal->Fill(totaldiff2);
            fDiffXvDiffY->Fill(xdiff_best,ydiff_best);
            fDiffDirCosX->Fill(lardirectionStart[0]-minostrack->ftrkdcosx);
            fDiffDirCosY->Fill(lardirectionStart[1]-minostrack->ftrkdcosy);
            fDiffDirCosZ->Fill(lardirectionStart[2]-minostrack->ftrkdcosz);
            fT962_DirCosX->Fill(lardirectionStart[0]);
            fT962_DirCosY->Fill(lardirectionStart[1]);
            fT962_DirCosZ->Fill(lardirectionStart[2]);
            fT962_StartZY->Fill(larStart[2],larStart[1]);
            fT962_EndZY->Fill(larEnd[2],larEnd[1]);
            fT962_StartXY->Fill(larStart[0],larStart[1]);
            fT962_EndXY->Fill(larEnd[0],larEnd[1]);
            fT962_StartZX->Fill(larStart[2],larStart[0]);
            fT962_EndZX->Fill(larEnd[2],larEnd[0]);
            fMinos_DirCosX->Fill(minostrack->ftrkdcosx);
            fMinos_DirCosY->Fill(minostrack->ftrkdcosy);
            fMinos_DirCosZ->Fill(minostrack->ftrkdcosz);
            fMinos_XY->Fill(100.0*minostrack->ftrkVtxX,100.0*minostrack->ftrkVtxY);
            if(minostrack->fcharge==1){
               fMinosErange_Pos->Fill(1000.0*minostrack->ftrkErange);
               fMinosMom_Pos->Fill(1000.0*minostrack->ftrkmom);
            }
            if(minostrack->fcharge==-1){
               fMinosErange_Neg->Fill(1000.0*minostrack->ftrkErange);
               fMinosMom_Neg->Fill(1000.0*minostrack->ftrkmom);
            }

            //make Association between T962 track and matched MINOS track
            util::CreateAssn(*this, evt, minostrack, lartrack, *assn);

         }//if(matchnumber!=999){

      }//loop over T962 tracks
                   
      evt.put(std::move(assn));//put Associations into event
        
      return;

   }
 
   //-------------------------------------------------
   bool MatchTracks::Compare(art::Ptr<recob::Track> lar_track, art::Ptr<t962::MINOS> minos_track,
                             double &xpred, double &ypred, double &rdiff, double &totaldiff, double &totaldiff2)
   {

      double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane 
                                           //(this minus number is the one we measured with Mitch)

      double x_offset=117.4; // previously 116.9;
      double y_offset=19.3; // previously  20.28;
      std::vector<double> larStart, larEnd;
      lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)
      //std::cout << "larStart = (" << larStart[0] << "," << larStart[1] << "," << larStart[2] << ")" << std::endl;
      //std::cout << "larEnd = (" << larEnd[0] << "," << larEnd[1] << "," << larEnd[2] << ")" << std::endl;

      double lardirectionStart[3];
      double lardirectionEnd[3];
      lar_track->Direction(lardirectionStart,lardirectionEnd);


      // double dz = D - larStart[2]+(100.0 * minos_track->ftrkVtxZ);//z-difference between end of T962 track and
      //begin of MINOS track...in centimeters

      double dz = D - larEnd[2]+(100.0 * minos_track->ftrkVtxZ);//z-difference between end of T962 track and
      //begin of MINOS track...in centimeters

      double l = dz/(lardirectionEnd[2]);//3-d distance between end of T962 track and begin of MINOS track

      double x_pred = l*lardirectionEnd[0]+larEnd[0];//predicted x-pos. of T962 track at z-position equal to
      //start of MINOS track
      double y_pred = l*lardirectionEnd[1]+larEnd[1];//predicted y-pos. of T962 track at z-position equal to
      //start of MINOS track

      double dx = 100.0*minos_track->ftrkVtxX - x_offset - x_pred;
      double dy = 100.0*minos_track->ftrkVtxY + y_offset - y_pred;
      
      if(100.0*minos_track->ftrkVtxZ > fdZ){
         //std::cout << "MINOS track #" << minos_track->ftrkIndex << " is too far into the detector.  Fail!" << std::endl;
         return false;//MINOS track starts too far into the detector 
      }

      if((x_pred + x_offset)>297.7 
         || (x_pred+x_offset)<-187.45 
         || (y_pred-y_offset)>181.71 
         || (y_pred-y_offset)<-100.29){
         return false;//predicted coordinates lie outside the MINOS detector boundaries
      }
      if(sqrt(dx*dx + dy*dy) > fdXY){
         return false;//predicted coordinates too far from XY-pos. at start of MINOS track
      }

      if( fabs(lardirectionEnd[0]-minos_track->ftrkdcosx) > fdCosx){
         return false;
      }
      if( fabs(lardirectionEnd[1]-minos_track->ftrkdcosy) > fdCosy){
         return false;
      }
      if( fabs(lardirectionEnd[2]-minos_track->ftrkdcosz) > fdCosz){
         return false;
      }

      xpred = dx;
      ypred = dy;
      rdiff = sqrt(dx*dx + dy*dy);
      //totaldiff is a measure of the agreement between the ArgoNeuT projected track and the MINOS track based on radial distance and angle. totaldiff= rdiff/cos(theta)  
      totaldiff=fabs(rdiff/((lardirectionEnd[0]*minos_track->ftrkdcosx)+(lardirectionEnd[1]*minos_track->ftrkdcosy)+(lardirectionEnd[2]*minos_track->ftrkdcosz)));

      if(totaldiff2>totaldiff){
         totaldiff2 = totaldiff;//new best match found
         return true;
      }
      else return false;//found match, but not best match
  
   }

   //--------------------------------------------------
   bool MatchTracks::EndsOnBoundary(art::Ptr<recob::Track> lar_track)
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

}//namespace match
