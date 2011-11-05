////////////////////////////////////////////////////////////////////////
/// \file  ShowerReco.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author:  biagio, ART port: echurch, detector agnostic: andrzejs
////////////////////////////////////////////////////////////////////////
#ifndef SHOWERRECO_H
#define SHOWERRECO_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework

#include <vector>
#include <string>

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TTree.h"

#include "RecoBase/Cluster.h"

namespace shwf {

  class ShowerReco : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerReco(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~ShowerReco();                               /**Destructor*/
    void beginJob();                                     
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Actual routine that reconstruct the shower*/
   
  //  int    Get3Daxis(float thetaI, float thetaC, float Wire_vertexI, float Wire_vertexC, float Time_vertex); // in rad

    int Get3DaxisN(int iplane0,int iplane1);

  //  int Get3Daxis_coords();

    //int    Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt); // in rad
   // void   GetVertex(art::Event& evt); // Get shower vertex
   // void   GetVertexN(art::Event& evt); // Get shower vertex

    void   GetVertexAndAnglesFromCluster(art::Ptr< recob::Cluster > clust,unsigned int plane); // Get shower vertex and slopes.

    void   GetPitchLength(); //Get pitch length of both planes
  //  void   AngularDistributionI(art::PtrVector < recob::Hit> hitlistInd); //Get angular distribution of the shower (Induction plane)
  //  void   AngularDistribution(art::PtrVector < recob::Hit> hitlist);  //Get angular distribution of the shower (Collection plane) 
    void   Get2DVariables(art::PtrVector < recob::Hit> hitlist);   

    void   FitAngularDistributions(int plane); 
 //   void   LongTransEnergyI(art::PtrVector < recob::Hit> hitlistInd); //Longtudinal and transverse enegry of the shower (Induction plane)
    void   LongTransEnergy(art::PtrVector < recob::Hit> hitlist); //Longtudinal and transverse enegry of the shower (Collection plane)
   double ProjectedLength(unsigned int plane ) 	 const;



    
//    float LastWire[2];  // last wire of the shower
//    float LastTime[2];  // last t_hit of the shower
    
    // 2D slope and intercept of the shower axis
    float slope[3];       // in cm, cm
    float intercept[3];   // in cm, cm
    float slope_wt[3];       // in wire, tick
    float intercept_wt[3];   // in wire,tick


      
    //double fTotCharge;  

   // float wireI1; // Wire Vertex position (Induction)
  //  float timeI1; // Time Vertex position (Induction)
  //  float wire1C; // Wire Vertex poistion (Collection)
   // float time1C; // Time Vertex position (Collection)
    const static float pi            = 3.141519;
    
    float ftimetick; // time sample in us
    const static float presamplings  = 60.;
    const static float fdriftvelocity =  0.157;
    const static int    alpha           = 5;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)

 private:

  double BirksCorrection(double dQdx_e);

   int fRun,fEvent,fSubRun;
   
   
//input labels:
  
  std::string fClusterModuleLabel;
  std::string fVertexCLusterModuleLabel;
    std::string fMCGeneratorLabel;
    std::string fLarGeantlabel;
  int fUseMCVertex;

 

  double xyz_vertex[3];


  float IdEdx4cm; // dedx of the first 4cm of the shower
  float CdEdx4cm; //dedx of the first 4cm of the shower 
  double CdEdx4cm_corr; //dedx of the first 4cm of the shower 
  double totCnrg,totCnrg_corr;
  double fMean_wire_pitch ;   // wire pitch in cm

//     std::vector<double> fOmega_Mean;    // Mean value of the 2D angular distribution (0=Ind - 1=Coll) cm,cm
//     std::vector<double> fOmega_RMS;     // RMS of the 2D angular distribution  (0=Ind - 1=Coll) cm, cm
// 
//     std::vector<double> fOmega_Mean_reb;    // Mean value of the 2D angular Rebinned by 4
//     std::vector<double> fOmega_RMS_reb;     // RMS of the 2D angular distribution  Rebinned by 4
//     std::vector<double> fOmega_Mean_Mean;    // Mean value of the 2D angular use mean instead of maximum
//  //   std::vector<double> fOmega_Mean_RMS;     // RMS of the 2D angular distribution use mean instead of maximum

   

    std::vector<double> fOmega_wt_Mean; // Mean value of the angular distribution (0=Ind - 1=Coll) wire,time
    std::vector<double> fOmega_wt_RMS;  // RMS of the angular distribution  (0=Ind - 1=Coll) wire,time

  std::vector<double> fTotChargeADC;   //Total charge in ADC/cm for each plane
  std::vector<double> fTotChargeMeV;  //Total charge in MeV/cm for each plane

  std::vector<double> fChargeADC_4cm;   //Initial charge in ADC/cm for each plane first 4cm
  std::vector<double> fChargeMeV_4cm;  //initial charge in MeV/cm for each plane first 4cm

  std::vector<double> fChargeADC_6cm;   //Initial charge in ADC/cm for each plane first 6cm
  std::vector<double> fChargeMeV_6cm;  //initial charge in MeV/cm for each plane first 6cm

  std::vector<double> fChargeADC_8cm;   //Initial charge in ADC/cm for each plane first 8cm
  std::vector<double> fChargeMeV_8cm;  //initial charge in MeV/cm for each plane first 8cm

  std::vector<double> fChargeADC_10cm;   //Initial charge in ADC/cm for each plane first 10cm
  std::vector<double> fChargeMeV_10cm;  //initial charge in MeV/cm for each plane first 10cm

  std::vector<std::vector<double> > fDistribChargeADC;  //vector with the first De/Dx points ADC
  std::vector<std::vector<double> > fDistribChargeMeV;  //vector with the first De/Dx points converted energy
  std::vector<std::vector<double> > fDistribChargeposition;  //vector with the first De/Dx points' positions 

std::vector<std::vector<double> > fSingleEvtAngle;  //vector with the first De/Dx points
std::vector<std::vector<double> > fSingleEvtAngleVal;  //vector with the first De/Dx points


    std::vector<unsigned int> fWire_vertex;  // wire coordinate of vertex for each plane
    std::vector<double> fTime_vertex;  // time coordinate of vertex for each plane

    std::vector<unsigned int> fWire_last;  // wire coordinate of last point for each plane
    std::vector<double> fTime_last;  // time coordinate of last point for each plane

 std::vector<unsigned int> fChannel_vertex;  // channel coordinate of vertex for each plane
 std::vector<unsigned int> fChannel_last;  // channel coordinate of vertex for each plane

   std::vector< double > fThetaN, fPhiN; // 3d angles calculated using a new "det independent" method
   std::vector< double > fThetaN_ang, fPhiN_ang; // 3d angles calculated using a new "det"// independent method degrees

//  std::vector< double > fThetaNC, fPhiNC; // 3d angles calculated using a new "det independent" method
//  std::vector< double > fThetaNC_ang, fPhiNC_ang; // 3d angles calculated using a new "det"
 //std::vector< double > fXvertex,fYvertex,fZvertex;
 //std::vector< double > fXlast,fYlast,fZlast;

       
  //  std::vector<double> fPitch;
    std::vector<double> fNPitch;


  
  unsigned int tpc;    //tpc type
  unsigned int fNPlanes; // number of planes  
  unsigned int fNAngles;
  
  
  TH1F* fh_theta[3]; 
  TH1F* fh_theta_wt[3]; 
  TH1F* fh_dedx[3]; 
  TH1F* fsh_nrg[3];
  TH1F* fsh_Tnrg[3];
  TH1F* fsh_long_hit[3];
  TTree* ftree_shwf;
TH1F *fh_phi_act;
TH1F *fh_thet_act;
  // Save enegry in a file
  //std::ofstream myfile;
  
  // Variables to get the angular distribution of the hits
 // float AI, BI, thetaI; // Induction  plane
  

 const static double calFactor=700.54;  // in ADC/fC
 const static double eCharge=1.6e-4;  // electron charge in fC


  }; // class ShowerReco

}

#endif // SHOWERRECO_H

