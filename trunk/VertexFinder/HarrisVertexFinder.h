////////////////////////////////////////////////////////////////////////
/// \file  HarrisVertexFinder.h
/// \brief Module to find vertices
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef HarrisVertexFinder_H
#define HarrisVertexFinder_H

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>

class TH2F;

///vertex reconstruction
namespace vertex {
   
 class HarrisVertexFinder :  public art::EDProducer {
    
  public:
    
    explicit HarrisVertexFinder(fhicl::ParameterSet const& pset); 
    virtual ~HarrisVertexFinder();        
    void beginJob();

    void produce(art::Event& evt);
    
  private:

    double Gaussian(int x, int y, double sigma);
    double GaussianDerivativeX(int x, int y);
    double GaussianDerivativeY(int x, int y);
  
    std::string fDBScanModuleLabel;
    int         fTimeBins;
    int         fMaxCorners;
    double      fGsigma;
    int         fWindow;
    double      fThreshold;
    int         fSaveVertexMap;
    TH2F *fNoVertices;

  };
    
}



#endif // HarrisVertexFinder_H
