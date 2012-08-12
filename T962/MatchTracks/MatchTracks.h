////////////////////////////////////////////////////////////////////////
//
// Compare ArgoNeuT/MINOS track collections and find any common tracks
//
// soderber@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef MATCHTRACKS_H
#define MATCHTRACKS_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"

#include "RecoBase/Track.h"
#include "T962/T962_Objects/MINOS.h"


#include <string>
#include <vector>


class TH1D;
class TH2D;

namespace match {
   
   class MatchTracks : public art::EDProducer {
      
   public:
      
      explicit MatchTracks(fhicl::ParameterSet const& pset); 
      virtual ~MatchTracks();      
      void beginJob();
      void produce(art::Event& evt);
      void reconfigure(fhicl::ParameterSet const& pset);
      
   private:

      TH1D* fDiffDirCosX;
      TH1D* fDiffDirCosY;
      TH1D* fDiffDirCosZ;
      TH1D* fDiffR;
      TH1D* fDiffTotal;
      TH2D* fDiffXvDiffY;
      TH1D* fT962_DirCosX;
      TH1D* fT962_DirCosY;
      TH1D* fT962_DirCosZ;
      TH2D* fT962_StartZY;
      TH2D* fT962_EndZY;
      TH2D* fT962_StartXY;
      TH2D* fT962_EndXY;
      TH2D* fT962_StartZX;
      TH2D* fT962_EndZX;
      TH1D* fT962_Ntracks;
      TH1D* fMINOS_Ntracks; 
      TH1D* fMinos_DirCosX;
      TH1D* fMinos_DirCosY;
      TH1D* fMinos_DirCosZ;
      TH2D* fMinos_XY;
      TH1D* fMinosErange_Pos;
      TH1D* fMinosErange_Neg;
      TH1D* fMinosMom_Pos;
      TH1D* fMinosMom_Neg;
      TH2D* fDiffXvD;
      TH2D* fDiffYvD;

      std::string fLarTracks_label;
      std::string fMinosTracks_label;
   
      double fdZ;    //max. z distance between projected T962 track and MINOS track
      double fdXY;   //max. radial distance between projected T962 track and MINOS track
      double fdCosx; //tolerance in dircosX matching
      double fdCosy; //tolerance in dircosY matching
      double fdCosz; //tolerance in dircosZ matching
      double fdBoundary; //distance from a boundary to be considered a track that "ends on a boundary"
      
   protected: 
      
      //method to compare track kinematics
      bool Compare(art::Ptr<recob::Track> lar_track, art::Ptr<t962::MINOS> minos_track, 
                   double &xpred, double &ypred, double &rdiff, double &totaldiff, double &totaldiff2); 
      bool EndsOnBoundary(art::Ptr<recob::Track> lar_track);

   }; // class MatchTracks
   
}

#endif // MATCHTRACKS_H
