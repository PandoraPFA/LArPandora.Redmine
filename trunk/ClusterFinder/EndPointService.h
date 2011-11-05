////////////////////////////////////////////////////////////////////////
/// \file  EndPointService.h
/// \brief Service to find 2D endpoints
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef EndPointService_H
#define EndPointService_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "TMath.h"
#include <vector>
#include <string>

namespace recob { 
  class Cluster;
  class EndPoint2D; 
}

namespace cluster {
   
  ///Service to find 2D end points
 class EndPointService {
    
  public:
    
    explicit EndPointService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg); 
    virtual ~EndPointService();        

    size_t EndPoint(art::PtrVector<recob::Cluster>& clusIn, std::vector<recob::EndPoint2D>& vtxcol);
    
  private:

    double Gaussian(int x, int y, double sigma);
    double GaussianDerivativeX(int x, int y);
    double GaussianDerivativeY(int x, int y);
    void VSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy);

    
    int          fTimeBins;
    int          fMaxCorners;
    double       fGsigma;
    int          fWindow;
    double       fThreshold;
    int          fSaveVertexMap;
  };
    
}



#endif // EndPointService_H
