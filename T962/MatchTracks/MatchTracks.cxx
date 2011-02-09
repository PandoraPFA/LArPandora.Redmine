////////////////////////////////////////////////////////////////////////
//
//   MatchTracks Class
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/OrphanHandle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "TFile.h"
#include "TH1D.h"


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
   MatchTracks::MatchTracks(fhicl::ParameterSet const& pset) : 
      fLarTracks_label(pset.get< std::string >("lartracks")),
      fMinosTracks_label(pset.get< std::string >("minostracks")),
      fdZ(pset.get< double>("dZ")),
      fdXY(pset.get< double>("dXY")),
      fdCosx(pset.get< double >("dCosx")),
      fdCosy(pset.get< double >("dCosy")),
      fdCosz(pset.get< double >("dCosz"))
   {

   }

   //-------------------------------------------------
   MatchTracks::~MatchTracks()
   {
   
   }

   //-------------------------------------------------
   void MatchTracks::beginJob()
   {
      //get access to the TFile service
      art::ServiceHandle<art::TFileService> tfs;

      fDiffDirCosX = tfs->make<TH1D>("fDiffDirCosX","Diff: DirCosX", 100,-2.0,2.0);
      fDiffDirCosY = tfs->make<TH1D>("fDiffDirCosY","Diff: DirCosY", 100,-2.0,2.0);
      fDiffDirCosZ = tfs->make<TH1D>("fDiffDirCosZ","Diff: DirCosZ", 100,-2.0,2.0);

   }

   //-------------------------------------------------
   void MatchTracks::analyze(const art::Event& evt)
   {

      art::Handle< std::vector<recob::Track> > LarTrackHandle;
      evt.getByLabel(fLarTracks_label,LarTrackHandle);

      art::Handle< std::vector<t962::MINOS> > MinosTrackHandle;
      evt.getByLabel(fMinosTracks_label,MinosTrackHandle);

      std::cout << "#T962 Tracks = " << LarTrackHandle->size() 
                << " #MINOS Tracks = " << MinosTrackHandle->size() << std::endl;

      for(unsigned int i=0; i<LarTrackHandle->size();++i){
         art::Ptr<recob::Track> lartrack(LarTrackHandle,i);
         
         for(unsigned int j=0; j<MinosTrackHandle->size();++j){
            art::Ptr<t962::MINOS> minostrack(MinosTrackHandle,j);

            bool match = Compare(lartrack,minostrack);
            if(match){
               std::cout << "Match!" << std::endl;
               double larStart[3];
               double larEnd[3];
               lartrack->Direction(larStart,larEnd);
               fDiffDirCosX->Fill(larStart[0]-minostrack->ftrkdcosx);
               fDiffDirCosY->Fill(larStart[1]-minostrack->ftrkdcosy);
               fDiffDirCosZ->Fill(larStart[2]-minostrack->ftrkdcosz);
            }
         }
      }
 
      return;
   }
 
   //-------------------------------------------------
   bool MatchTracks::Compare(art::Ptr<recob::Track> lar_track, art::Ptr<t962::MINOS> minos_track)
   {

       double D=(90*0.5)+(42.4*2.54)-5.588; //distance from the front (upstream) of the TPC to the 1st Minos plane 
                                           //(this minus number is the one we measured with Mitch)

      double x_offset=116.9; // previously 118;
      double y_offset=20.28; //previously  19;

      std::vector<double> larStart, larEnd;
      lar_track->Extent(larStart,larEnd);//put xyz coordinates at begin/end of track into vectors(?)


      double lardirectionStart[3];
      double lardirectionEnd[3];
      lar_track->Direction(lardirectionStart,lardirectionEnd);


      double dz = D - larStart[2]+(100.0 * minos_track->ftrkVtxZ);//z-difference between end of T962 track and
                                                                  //begin of MINOS track...in centimeters

      double l = dz/(lardirectionEnd[2]);//3-d distance between end of T962 track and begin of MINOS track

      double x_pred = l*lardirectionEnd[0]+larEnd[0];//predicted x-pos. of T962 track at z-position equal to
                                                       //start of MINOS track
      double y_pred = l*lardirectionEnd[1]+larEnd[1];//predicted y-pos. of T962 track at z-position equal to
                                                       //start of MINOS track

      double dx = 100.0*minos_track->ftrkVtxX - x_pred;
      double dy = 100.0*minos_track->ftrkVtxY - y_pred;
      
      
      if(100.0*minos_track->ftrkVtxZ > fdZ) return false;//MINOS track starts too far into the detector 
      
      if((x_pred + x_offset)>297.7 
         || (x_pred+x_offset)<-187.45 
         || (y_pred+y_offset)>181.71 
         || (y_pred+y_offset)<-100.29) return false;//predicted coordinates lie outside the MINOS detector boundaries

      if(sqrt(dx*dx + dy*dy) > fdXY) return false;//predicted coordinates too far from XY-pos. at start of MINOS track

      if( fabs(lardirectionStart[0]-minos_track->ftrkdcosx) > fdCosx) return false;
      if( fabs(lardirectionStart[1]-minos_track->ftrkdcosy) > fdCosy) return false;
      if( fabs(lardirectionStart[2]-minos_track->ftrkdcosz) > fdCosz) return false;
     
      return true;  
  
   }

}//namespace match
