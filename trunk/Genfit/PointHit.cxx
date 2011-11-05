// This Class' Header ------------------
#include "PointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "GeaneTrackRep2.h"
#include "RKTrackRep.h"
#include "SlTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

//ClassImp(PointHit)


genf::PointHit::~PointHit()
{}

genf::PointHit::PointHit()
  : SpacepointRecoHit(NparHitRep)
{}

genf::PointHit::PointHit(TVector3 point,double res)
  : SpacepointRecoHit(NparHitRep){

  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;
  fHitCov[2][2] = res*res;
  GFDetPlane d;


  fHitCoord[0][0] = gRandom->Gaus(point.X(),res);
  fHitCoord[1][0] = gRandom->Gaus(point.Y(),res);
  fHitCoord[2][0] = gRandom->Gaus(point.Z(),res);



}
genf::GFAbsRecoHit* genf::PointHit::clone(){
  return new PointHit(*this);
}


TMatrixT<Double_t>
genf::PointHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const genf::GeaneTrackRep2*>(stateVector) != NULL) ||
      (dynamic_cast<const genf::RKTrackRep*>(stateVector) != NULL)){
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);


    //   TMatrixT<Double_t> HMatrix(2,6);
    // WTF? Making this 2,5 not 2,6. EC, 3-Jan-2011.
    TMatrixT<Double_t> HMatrix(2,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    //HMatrix[0][5] = 0.;

    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;
    //HMatrix[1][5] = 0.;

    return HMatrix;
  }
  else if(dynamic_cast<const genf::SlTrackRep*>(stateVector)){
    TMatrixT<Double_t> HMatrix(2,4);
    
    HMatrix[0][0] = 1.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;

    
    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 1.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    return HMatrix;
  }

 else {
   std::cerr << "PointHit can only handle state"
	     << " vectors of type GeaneTrackRep2 -> abort" << std::endl;
   throw;
 }
 
}


