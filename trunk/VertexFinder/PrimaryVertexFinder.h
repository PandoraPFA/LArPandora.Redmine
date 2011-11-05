////////////////////////////////////////////////////////////////////////
/// \file  PrimaryVertexFinder.h
/// \brief Module to find vertices
///
/// \saima@ksu.edu
////////////////////////////////////////////////////////////////////////

#ifndef PrimaryVertexFinder_H
#define PrimaryVertexFinder_H

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include "RecoBase/recobase.h"

class TH1F;
class TH2F;

///vertex reconstruction
namespace vertex {
   
 class PrimaryVertexFinder :  public art::EDProducer {
    
  public:
    
    explicit PrimaryVertexFinder(fhicl::ParameterSet const& pset); 
    virtual ~PrimaryVertexFinder();        
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);

  private:

  
    std::string fTrackModuleLabel;
    double      fVertexWindow;
    double      StartPointSeperation(recob::SpacePoint sp1, recob::SpacePoint sp2);
    bool        IsInVertexCollection(int a, std::vector<std::vector<int> > vertex_collection);
    int         IndexInVertexCollection(int a, int b, std::vector<std::vector<int> > vertex_collection);
    bool        IsInNewVertex(int a, std::vector<int> newvertex);
    double      gammavalue(TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2);
    double      alphavalue(double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2);
    double      MinDist(double alpha, double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2);
    TVector3    PointOnExtendedTrack(double alphagamma, TVector3 startpoint,  TVector3 dircos);
    TH2F*       fNoTracks;
    TH1F*       fLength_1stTrack;
    TH1F*       fLength_2ndTrack;
    TH1F*       fLength_3rdTrack;
    TH1F*       fLength_4thTrack;
    TH1F*       fLength_5thTrack;

  };
    
}



#endif // PrimaryVertexFinder_H
