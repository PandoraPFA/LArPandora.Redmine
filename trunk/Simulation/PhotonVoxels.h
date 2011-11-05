#include "TVector3.h"
#include "iostream"

#ifndef PhotonVoxels_h 
#define PhotonVoxels_h 1


namespace sim {


  class PhotonVoxel
  {
  public:
    PhotonVoxel(Double_t, Double_t, Double_t, Double_t, Double_t,\
		Double_t, int N = 0) ;
    PhotonVoxel();
    ~PhotonVoxel();

    TVector3 GetLowerCorner() const;
    TVector3 GetUpperCorner() const;
    TVector3 GetCenter() const;


  private:
    double xVoxelMin, xVoxelMax;
    double yVoxelMin, yVoxelMax;
    double zVoxelMin, zVoxelMax;

    int NPhotons;

  };


  class PhotonVoxelDef
  {
  public:
    PhotonVoxelDef(Double_t, Double_t, int, Double_t, Double_t,\
		   int, Double_t, Double_t, int );
    PhotonVoxelDef();
    ~PhotonVoxelDef();

    TVector3 GetRegionUpperCorner() const;
    TVector3 GetRegionLowerCorner() const;
    TVector3 GetSteps() const;

    TVector3 GetVoxelSize()const;

    int GetNVoxels() const;

    int GetVoxelID(TVector3) const;
    bool IsLegalVoxelID(int) const;

    PhotonVoxel GetPhotonVoxel(int ID) const;
    PhotonVoxel GetContainingVoxel(TVector3) const;

    bool operator==(const PhotonVoxelDef &rhs) const;
    bool operator!=(const PhotonVoxelDef &rhs) const 
      { return ! ((*this)==rhs); }


  private:
    TVector3 fLowerCorner;
    TVector3 fUpperCorner;
    int fxSteps, fySteps, fzSteps;

  };
}

#endif
