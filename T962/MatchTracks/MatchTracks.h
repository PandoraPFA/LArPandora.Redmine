////////////////////////////////////////////////////////////////////////
//
// Compare ArgoNeuT/MINOS track collections and find any common tracks
//
// soderber@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef MATCHTRACKS_H
#define MATCHTRACKS_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Persistency/Common/Ptr.h"

#include "RecoBase/Track.h"
#include "T962/T962_Objects/MINOS.h"

#include <string>
#include <vector>


class TH1D;

namespace match {
   
   class MatchTracks : public art::EDAnalyzer {
      
   public:
      
      explicit MatchTracks(fhicl::ParameterSet const& pset); 
      virtual ~MatchTracks();
      
      void beginJob();
      void analyze(const art::Event& evt);
      
   private:

      TH1D* fDiffDirCosX;
      TH1D* fDiffDirCosY;
      TH1D* fDiffDirCosZ;

      std::string fLarTracks_label;
      std::string fMinosTracks_label;
   
      double fdZ;    //max. z distance between projected T962 track and MINOS track
      double fdXY;   //max. radial distance between projected T962 track and MINOS track
      double fdCosx; //tolerance in dircosX matching
      double fdCosy; //tolerance in dircosY matching
      double fdCosz; //tolerance in dircosZ matching
      
   protected: 
      
      bool Compare(art::Ptr<recob::Track> lar_track, art::Ptr<t962::MINOS> minos_track); //method to compare track kinematics
      
   }; // class MatchTracks
   
}

#endif // MATCHTRACKS_H
