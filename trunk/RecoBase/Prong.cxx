////////////////////////////////////////////////////////////////////////////
// \version $Id: Prong.cxx,v 1.4 2010/06/10 16:21:31 antonm Exp $
//
// \brief Definition of prong object for LArSoft
//
// \author brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoBase/Prong.h"
#include "TMath.h"
#include "TVector3.h"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>

namespace recob{

  //----------------------------------------------------------------------
  Prong::Prong()
  {
  }

  //----------------------------------------------------------------------
  Prong::Prong(art::PtrVector<recob::Cluster> &clusters, 
	       std::vector<recob::SpacePoint> spacepoints) :
    fClusters(clusters),
    fSpacePoints(spacepoints)
  {
    mf::LogWarning("Prong") << "need to define fMaxTransverseWidth and fDistanceMaxWidth in ctor";
  }

  //----------------------------------------------------------------------
  Prong::~Prong()
  {
  }

  //----------------------------------------------------------------------
  // in the case of view == geo::kUnknown use all hits in all clusters
  art::PtrVector<recob::Hit> Prong::Hits(geo::View_t view) const
  {
    art::PtrVector<recob::Hit> hits;
    for(size_t i = 0; i < fClusters.size(); ++i){
      if(fClusters[i]->View() == view || view == geo::kUnknown){
	art::PtrVector<recob::Hit> hs = fClusters[i]->Hits();
        for(size_t ch = 0; ch < hs.size(); ++ch) hits.push_back(hs[ch]);
      }
    }

    return hits;
  }

  //----------------------------------------------------------------------
  art::PtrVector<recob::Cluster> Prong::Clusters(geo::View_t view) const
  {
    art::PtrVector<recob::Cluster> clusts;
    for(size_t i = 0; i < fClusters.size(); ++i){
      if(fClusters[i]->View() == view) clusts.push_back(fClusters[i]);
    }

    return clusts;
  }

  //----------------------------------------------------------------------
  void Prong::SetDirection(double *dcosStart, double *dcosEnd)
  {
    fDCosStart[0] = dcosStart[0];
    fDCosStart[1] = dcosStart[1];
    fDCosStart[2] = dcosStart[2];

    fDCosEnd[0] = dcosEnd[0];
    fDCosEnd[1] = dcosEnd[1];
    fDCosEnd[2] = dcosEnd[2];

    return;
  }

  //----------------------------------------------------------------------
  void Prong::Direction(double *dcosStart, double *dcosEnd) const
  {
    dcosStart[0] = fDCosStart[0];
    dcosStart[1] = fDCosStart[1];
    dcosStart[2] = fDCosStart[2];

    dcosEnd[0] = fDCosEnd[0];
    dcosEnd[1] = fDCosEnd[1];
    dcosEnd[2] = fDCosEnd[2];

    return;
  }

  //----------------------------------------------------------------------
  double Prong::ProjectedLength(geo::View_t view) const
  {
    double pitch = -1.;
    if(view == geo::kUnknown || view == geo::k3D){
      mf::LogWarning("RecoBaseProng")<< "Warning Prong::ProjectedLength :  no Pitch foreseen for view "<<view;
      return pitch;
    }
    else{
      art::ServiceHandle<geo::Geometry> geo;
      for(unsigned int i = 0; i < geo->Nplanes(); ++i){
	if(geo->Plane(i).View() == view){
	  double wirePitch = geo->WirePitch(0,1,i);
	  double angleToVert = geo->Plane(i).Wire(0).ThetaZ(false) - 0.5*TMath::Pi();
	  //(sin(angleToVert),cos(angleToVert)) is the direction perpendicular to wire
	  //(fDCosStart[1],fDDosStart[2]) is the direction of prong in the y-z plane
	  double cosgamma = TMath::Abs(TMath::Sin(angleToVert)*fDCosStart[1]+TMath::Cos(angleToVert)*fDCosStart[2]);
	  if (cosgamma>0) pitch = wirePitch/cosgamma;	  
	} // end if the correct view
      } // end loop over planes
    } // end if a reasonable view

    return pitch;
  }

  //----------------------------------------------------------------------
  void Prong::Extent(std::vector<double> &xyzStart, std::vector<double> &xyzEnd) const
  {
    // get the first and last hits
    double *Start = (double *)(fSpacePoints.front().XYZ());
    double *End   = (double *)(fSpacePoints.back().XYZ());

    xyzStart.clear();
    xyzEnd.clear();
    for(unsigned int i = 0;i<3;i++){
      xyzStart.push_back(*(Start+i));
      xyzEnd.push_back(*(End+i));
    }
    
    return;
  }

  //----------------------------------------------------------------------
  double Prong::MaxTransverseWidth(geo::Coord_t coordinate) const
  {
    if(coordinate      == geo::kX) return fMaxTransverseWidth[0];
    else if(coordinate == geo::kY) return fMaxTransverseWidth[1];

    mf::LogWarning("RecoBaseProng") << "supplied coordinate to MaxTransverseWidth is unknown, return X vaule";
    return fMaxTransverseWidth[0];
  }

    
  //----------------------------------------------------------------------
  double Prong::Theta() const
  {
      TVector3 vec(fDCosStart[0],fDCosStart[1],fDCosStart[2]);
      return vec.Theta();
  }

  //----------------------------------------------------------------------
  double Prong::Phi() const
  {
      TVector3 vec(fDCosStart[0],fDCosStart[1],fDCosStart[2]);
      return vec.Phi();
  }
    
    
  //----------------------------------------------------------------------
  // Print function...to be called by << operator.  Facilitates overriding
  // by inheriting Track/Shower classes.
  //
   std::ostream& Prong::Print(std::ostream& stream) const
   {
      double dcoss[3];
      double dcose[3];
      this->Direction(dcoss,dcose);
      stream << std::setprecision(5);
      stream << "Prong ID " << std::setw(5)  << std::right << this->ID() 
              << " #SpacePoints = " << std::setw(5)  << std::right << this->SpacePoints().size()
              << " #Clusters = " << std::setw(5) << std::right <<  this->fClusters.size();

      return stream;
     
   }

  //----------------------------------------------------------------------
  // ostream operator.  
  //
  std::ostream& operator<< (std::ostream& o, const Prong & a)
  {
     return a.Print(o);
  }


  //----------------------------------------------------------------------
  // < operator.  
  //
  bool operator < (const Prong & a, const Prong & b)
  {
    if(a.ID() != b. ID())
      return a.ID()<b.ID();

    return false; //They are equal
  }


}
