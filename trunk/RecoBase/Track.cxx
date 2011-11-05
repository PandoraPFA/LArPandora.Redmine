////////////////////////////////////////////////////////////////////////////
// \version $Id: Track.cxx,v 1.5 2010/02/15 20:32:46 brebel Exp $
//
// \brief Definition of track object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////

#include "RecoBase/Track.h"
#include <string>
#include <iostream>

namespace recob{

   //----------------------------------------------------------------------
   Track::Track() :
      Prong()
   {
   }
   
   //----------------------------------------------------------------------
   Track::Track(art::PtrVector<recob::Cluster> &clusters, std::vector<recob::SpacePoint> spacepoints) :
      Prong(clusters,spacepoints)
   {
   }
   //----------------------------------------------------------------------
   Track::~Track()
   {
   }
   
   //----------------------------------------------------------------------
   // ostream operator.  
   //
   std::ostream& Track::Print(std::ostream& stream) const
   {
      double dcoss[3];
      double dcose[3];
      this->Direction(dcoss,dcose);
      stream << std::setiosflags(std::ios::fixed) << std::setprecision(3);
      stream << "Track ID " << std::setw(4)  << std::right << this->ID() 
             << " : #SpacePoints = " << std::setw(5)  << std::right << this->SpacePoints().size()
             << " #Clusters = " << std::setw(3) << std::right << this->fClusters.size() 
             << " StartCosines : ( " ;
      for(int i = 0; i<3; ++i) stream << std::setw(6) << std::right << dcoss[i] << " ";
      stream << ")  EndCosines : ( ";
      for(int i = 0; i<3; ++i) stream << std::setw(6) << std::right << dcose[i] << " ";
      stream << ")" ;
      stream << " Theta = " << std::setw(6) << std::right << this->Theta() << " ";
      stream << " Phi = " << std::setw(6) << std::right << this->Phi() ;
      return stream;

   }
   

}
