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
      fdcosx(pset.get< double >("dcosx")),
      fdcosy(pset.get< double >("dcosy")),
      fdcosz(pset.get< double >("dcosz"))
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
      double lardirectionStart[3];
      double lardirectionEnd[3];
      lar_track->Direction(lardirectionStart,lardirectionEnd);
      if( fabs(lardirectionStart[0]-minos_track->ftrkdcosx) > fdcosx) return false;
      if( fabs(lardirectionStart[1]-minos_track->ftrkdcosy) > fdcosy) return false;
      if( fabs(lardirectionStart[2]-minos_track->ftrkdcosz) > fdcosz) return false;
     
      return true;  
  
   }

}//namespace match
