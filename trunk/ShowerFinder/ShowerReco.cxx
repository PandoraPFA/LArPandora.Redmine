////////////////////////////////////////////////////////////////////////
//
// \file ShowerReco.cxx
//
// biagio.rossi@lhep.unibe.ch   (FWMK : argoneut specific)
// thomas.strauss@lhep.unibe.ch (ART  : general detector)
//
// andrzej.szelc@yale.edu (port to detector agnostic version)
//
// This algorithm is designed to reconstruct showers
// 
///////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

// LArSoft includes

#include "ShowerFinder/ShowerReco.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"



#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/summary.h"


// ***************** //


shwf::ShowerReco::ShowerReco(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Shower> >();
}

void shwf::ShowerReco::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel = pset.get< std::string >("ClusterModuleLabel");
   //fHoughLineModuleLabel=pset.get<std::string > ("HoughLineModuleLabel");
  fVertexCLusterModuleLabel=pset.get<std::string > ("VertexClusterModuleLabel");
 // fMCGeneratorLabel=pset.get<std::string > ("MCGeneratorLabel");
  //fLarGeantlabel=pset.get<std::string > ("LarGeantlabel");     
  //fUseMCVertex=pset.get<int > ("UseMCVertex");
 

}

// ***************** //
shwf::ShowerReco::~ShowerReco()
{
}

namespace shwf {
struct SortByWire 
{
  bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
  { return 
      h1.Wire()->RawDigit()->Channel() < 
      h2.Wire()->RawDigit()->Channel() ;
  }
};
}

// ***************** //
void shwf::ShowerReco::beginJob()
{



  /** Get Geometry*/
  art::ServiceHandle<geo::Geometry> geo;
  fNPlanes = geo->Nplanes();
  fMean_wire_pitch = geo->WirePitch(0,1,0);    //wire pitch in cm

  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  art::ServiceHandle<util::DetectorProperties> detp;
  ftimetick=detp->SamplingRate()/1000.;
  
  std::cout << "------------ timetick? " << ftimetick << std::endl << std::endl << std::endl;

  
  
  /** Create Histos names*/
  char tit_dedx[128] = {0};
  char tit_h_theta[128] = {0};
  char tit_h_phi_act[128] = {0};
  char tit_h_thet_act[128] = {0};
  char sh_long_tit[128] = {0};
  char sh_tit[128] = {0};
  char shT_tit[128] = {0};
  
  int nbins;
  
   /**Histos for the actual distribution of the angle in 3D*/
     
    sprintf(tit_h_phi_act,"fh_phi_act");
    fh_phi_act = tfs->make<TH1F>(tit_h_phi_act,"3D phi distribution",720,-180., 180.);

    sprintf(tit_h_thet_act,"fh_thet_act");
    fh_thet_act = tfs->make<TH1F>(tit_h_thet_act,"3D Theta distribution",720,-180., 180.);

  //int tt=0;
  
  
  for(unsigned int i=0;i<fNPlanes;++i){


    //    sprintf(&tit_dedx[0],"fh_dedx_%.4i_%.4i_%i",i);
    sprintf(&tit_dedx[0],"fh_dedx_._%i",i);
    fh_dedx[i] = tfs->make<TH1F>(tit_dedx,"dEdx vs distance from vertex",120,0, 40);

    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_theta_%i",i);
    fh_theta[i] = tfs->make<TH1F>(tit_h_theta,"Theta distribution",720,-180., 180.);

    
      /**Histos for the longitudinal energy distribution of the shower */
    sprintf(&sh_tit[0],"fsh_nrg1_%i",i);                   /**number of wires used,min wire,max_wire you need for the anlysis*/
    fsh_nrg[i] = tfs->make<TH1F>(sh_tit,"energy reco",240,0.,240*fMean_wire_pitch);

    /**Histos for the transverse energy distribution of the shower*/
    sprintf(&shT_tit[0],"fshT_nrg1_%i",i);                   /**units are ticks most lickely, but how did Biagio get size of it???*/
    fsh_Tnrg[i] = tfs->make<TH1F>(shT_tit,"energy reco",80,-40.,40.);
  
    /**Histos for the Transverse HIT distribution of the shower*/
    nbins = (int)(240*fMean_wire_pitch);
    sprintf(&sh_long_tit[0],"fsh_long_hit_%i",i);                           /**nbins,min wire,max_wire you need for the analysis*/
    fsh_long_hit[i] = tfs->make<TH1F>(sh_long_tit,"longitudinal hit reco",nbins, 0.,     240*fMean_wire_pitch);
  }

  ftree_shwf =tfs->make<TTree>("ShowerReco","Results");/**All-knowing tree with reconstruction information*/
  
  // ftree_shwf->Branch("ftheta_Mean","std::vector<double>",&fTheta_Mean);
  // ftree_shwf->Branch("ftheta_RMS","std::vector<double>",&fTheta_RMS );
   
   ftree_shwf->Branch("run",&fRun,"run/I");
    ftree_shwf->Branch("subrun",&fSubRun,"subrun/I");
   ftree_shwf->Branch("event",&fEvent,"event/I");
   ftree_shwf->Branch("nplanes",&fNPlanes,"nplanes/I");
    ftree_shwf->Branch("nangles",&fNAngles,"nangles/I");
   
   ftree_shwf->Branch("fthetaN","std::vector<double>",&fThetaN_ang);
   ftree_shwf->Branch("fphiN","std::vector<double>",&fPhiN_ang);

//    ftree_shwf->Branch("fthetaNC","std::vector<double>",&fThetaNC_ang);
//    ftree_shwf->Branch("fphiNC","std::vector<double>",&fPhiNC_ang);
   
   
   ftree_shwf->Branch("ftotChargeADC","std::vector<double>",&fTotChargeADC);
   ftree_shwf->Branch("ftotChargeMeV","std::vector<double>",&fTotChargeMeV);
   

   ftree_shwf->Branch("NPitch","std::vector<double>", &fNPitch);
 //  ftree_shwf->Branch("Pitch","std::vector<double>", &fPitch);

// this should be temporary - until the omega is sorted out.
  
   
   ftree_shwf->Branch("ChargeADC_4cm","std::vector<double>",&fChargeADC_4cm);
   ftree_shwf->Branch("ChargeMeV_4cm","std::vector<double>",&fChargeMeV_4cm);
   ftree_shwf->Branch("ChargeADC_6cm","std::vector<double>",&fChargeADC_6cm);
   ftree_shwf->Branch("ChargeMeV_6cm","std::vector<double>",&fChargeMeV_6cm);
   ftree_shwf->Branch("ChargeADC_8cm","std::vector<double>",&fChargeADC_8cm);
   ftree_shwf->Branch("ChargeMeV_8cm","std::vector<double>",&fChargeMeV_8cm);
   ftree_shwf->Branch("ChargeADC_10cm","std::vector<double>",&fChargeADC_10cm);
   ftree_shwf->Branch("ChargeMeV_10cm","std::vector<double>",&fChargeMeV_10cm);


   ftree_shwf->Branch("ChargedistributionADC","std::vector<std::vector<double>>",&fDistribChargeADC);
   
   ftree_shwf->Branch("ChargedistributionMeV","std::vector<std::vector<double>>",&fDistribChargeMeV);
   
   ftree_shwf->Branch("ChargedistributionPosition","std::vector<std::vector<double>>",&fDistribChargeposition);
 
 
}

// ***************** //
void shwf::ShowerReco::produce(art::Event& evt)
{ 

  /**TODO THIS VALUES SHOULD BE PARAMETERS OF THE MODULE, evtl from database
  *double Efield_SI        =  0.7;     // Electric Field between Shield and Induction planes in kV/cm
  *double Efield_IC        =  0.9;     // Electric Field between Induction and Collection planes in kV/cm
  *double Temperature      = 87.6;  // LAr Temperature in K
  *check if there can be a replacement later for the product, not needed now
  *double timepitch        = fdriftvelocity*ftimetick;                  //time sample (cm) */

  /** Get Geometry */
  art::ServiceHandle<geo::Geometry> geo;
  //art::ServiceHandle<util::LArProperties> larprop;
  fNPlanes = geo->Nplanes();
  //fdriftvelocity = larprob->DriftVelocity(Efield_SI,Temperature);
  
  //calculate factorial for number of angles
  int fact=1;
  for (unsigned int i=1; i<=fNPlanes; i++)
    fact*=i;
  
  fNAngles=fact/2;

  fPhiN.resize(0);
  fThetaN.resize(0);
 
  fPhiN_ang.resize(0);
  fThetaN_ang.resize(0);
 
  
fDistribChargeADC.resize(fNPlanes); 
fDistribChargeMeV.resize(fNPlanes);
fDistribChargeposition.resize(fNPlanes);

 for(unsigned int ii=0;ii<fNPlanes;ii++)
  {  fDistribChargeADC[ii].resize(0);  //vector with the first De/Dx points
  fDistribChargeMeV[ii].resize(0);  //vector with the first De/Dx points
  fDistribChargeposition[ii].resize(0);  //vector with the first De/Dx points' positions 
  }


      // fPitch.resize(fNPlanes); 
    	fNPitch.resize(fNPlanes,-1); 
	fWire_vertex.resize(fNPlanes,-1);
	fTime_vertex.resize(fNPlanes,-1);
        fWire_last.resize(fNPlanes,-1);
	fTime_last.resize(fNPlanes,-1);
 	fTotChargeADC.resize(fNPlanes,0); 
	fTotChargeMeV.resize(fNPlanes,0);  
      //  fChannel_vertex.resize(fNPlanes,-1);
      //  fChannel_last.resize(fNPlanes,-1);

fChargeADC_4cm.resize(fNPlanes,0);   //Initial charge in ADC/cm for each plane angle calculation 4cm
  fChargeMeV_4cm.resize(fNPlanes,0);  //initial charge in MeV/cm for each angle calculation first 4cm
  fChargeADC_6cm.resize(fNPlanes,0);   //Initial charge in ADC/cm for each angle calculation first 6cm
  fChargeMeV_6cm.resize(fNPlanes,0);  //initial charge in MeV/cm for each angle calculation first 6cm
  fChargeADC_8cm.resize(fNPlanes,0);   //Initial charge in ADC/cm for each angle calculationfirst 8cm
  fChargeMeV_8cm.resize(fNPlanes,0);  //initial charge in MeV/cm for each angle calculation first 8cm
  fChargeADC_10cm.resize(fNPlanes,0);   //Initial charge in ADC/cm for each angle calculation first 10cm
  fChargeMeV_10cm.resize(fNPlanes,0);  //initial charge in MeV/cm for each angle calculation first 10cm

//       


  /**Get Clusters*/
  std::cout << "************ What I'm getting out " << fClusterModuleLabel << " " << std::endl;
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);


   std::vector< art::PtrVector < recob::Hit> > hitlist_all;
   //art::PtrVector < recob::Hit> hitlistInd;
   hitlist_all.resize(fNPlanes);

 
    //std::auto_ptr<std::vector<recob::Shower> > Shower3DVector(new std::vector<recob::Shower>);

  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {

      art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      
      //get vertex position and slope information to start with:
     
      art::PtrVector<recob::Hit> hitlist;
      hitlist = cl->Hits();
      hitlist.sort(shwf::SortByWire());
      unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire


      
      geo->ChannelToWire((*hitlist.begin())->Wire()->RawDigit()->Channel(),t,p,w);
      
      GetVertexAndAnglesFromCluster( cl,p);

      
      
      for(art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin(); a != hitlist.end();  a++) //loop over cluster hits
      {
	c=(*a)->Wire()->RawDigit()->Channel(); 
	geo->ChannelToWire(c,t,p,w);

        //geo->Plane(i).View()

	//std::cout << "+++++++ planeview in hit list   " <<  geo->Plane(p,t).View() << std::endl; 

          hitlist_all[p].push_back(*a);

      }

    } // End loop on clusters.

 fRun = evt.id().run();
 fSubRun = evt.id().subRun();
 fEvent = evt.id().event();



   for(unsigned int i=0;i<fNPlanes;i++)
      {
       hitlist_all[i].sort(shwf::SortByWire());
      // Get2DVariables(hitlist_all[i]);
       }
 
 
  
///////////////// !!! this has to be plane independent
 // Get3Daxis(fOmega_Mean[1], fOmega_Mean[0], fWire_vertex[1], fWire_vertex[0], fTime_vertex[0]);
 
 for(unsigned int ii=0;ii<fNPlanes-1;ii++)
    for(unsigned int ij=ii+1;ij<fNPlanes;ij++)
      Get3DaxisN(ii,ij);
  
    
   /* Get3DaxisN(0,1);
if(fNPlanes>2)
  {	
  Get3DaxisN(0,2);
  Get3DaxisN(1,2);
  }*/	
  
// Get3Daxis_coords();  


for(unsigned int i=0;i<fNPlanes;i++)
  LongTransEnergy(hitlist_all[i]); //Longitudinal and Transverse energy profile of the Shower induction
  
// LongTransEnergy(hitlistCol); //Longitudinal and Transverse energy profile of the Shower induction
 ///////////////// !!! this has to be plane independent 
 

//////create spacepoints, and direction cosines for Shower creation

//std::vector< recob::SpacePoint > 	spacepoints = std::vector<recob::SpacePoint>()	


         
// make an art::PtrVector of the clusters
 art::PtrVector<recob::Cluster> prodvec;
 for(unsigned int i = 0; i < clusterListHandle->size(); ++i){
   art::Ptr<recob::Cluster> prod(clusterListHandle, i);
   prodvec.push_back(prod);
 }

//create a singleSpacePoint at vertex.
std::vector< recob::SpacePoint > spcpts;

//art::PtrVector< recob::Hit > hits;

//for(int pl=0;pl<fNPlanes;pl++)
 //{
 //for each plane create or find hit closest to the vertex... 
 //recob::Hit  sinhit(art::Ptr< recob::Wire > &wire, 0., 0., 0., 0., fTime_vertex[pl], 0., 0., 0., 0., 0., 0., 0.)

 //}

//double xyz[3]={};

recob::SpacePoint singlepoint;
singlepoint.SetXYZ(xyz_vertex);
spcpts.push_back(singlepoint);


recob::Shower  singShower(prodvec,spcpts);


// TBD determine which angle to use for the actual shower
double fPhi=fPhiN_ang[0];
double fTheta=fThetaN_ang[0];

double dcosstart[3]={TMath::Cos(fPhi*pi/180)*TMath::Sin(fTheta*pi/180),TMath::Cos(fTheta*pi/180),TMath::Sin(fPhi*pi/180)*TMath::Sin(fTheta*pi/180)};


singShower.SetDirection(dcosstart,dcosstart);
singShower.SetID(1);

std::auto_ptr<std::vector<recob::Shower> > Shower3DVector(new std::vector<recob::Shower>);
Shower3DVector->push_back(singShower);

//get direction cosines and set them for the shower




//


  /**Fill the output tree with all information */
    ftree_shwf->Fill();

//for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane)
  //fh_theta[iplane]->Write(Form("fh_theta_%d_%d",iplane,evt.id().event()));
  // This needs work, clearly.  
  //for(int p=0;p<2;p++)Shower3DVector->push_back(shower);
  evt.put(Shower3DVector);

}



// void shwf::ShowerReco::Get2DVariables(art::PtrVector < recob::Hit> hitlist) {  
// 
//   art::ServiceHandle<geo::Geometry> geom;
//   // only needed for drawing the axis of the shower in the event display
//   
//   unsigned int channel;
//   double  wire_cm, time_cm;
// 
// 
//  double AC, BC, omega; 
//   double time;
//   unsigned int wire,plane;
// 
// for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
//     art::Ptr<recob::Hit> theHit = (*hitIter);
//     time = theHit->PeakTime() ;  
//     //time_C -= (presamplings+10.1);
//     art::Ptr<recob::Wire> theWire = theHit->Wire();
//     channel = theWire->RawDigit()->Channel();
//     geom->ChannelToWire(channel, tpc, plane, wire);
//     //    if(time_C<1020)continue;
//     wire_cm = wire * fMean_wire_pitch; //in cm
//     time_cm = time *ftimetick*fdriftvelocity; //in cm
//    
//  /*  if(hitIter == hitlist.begin())
// 	{ fWire_vertex[plane] = fWire_vertex[plane] ; //in cm
//           fTime_vertex[plane] = fTime_vertex[plane] ; //in cmm        
//         }  */  
// 
//    
//     // moving to polar coordinates
//     BC = (wire_cm - fWire_vertex[plane]* fMean_wire_pitch)+fMean_wire_pitch; //in cm
//     AC = (time_cm - fTime_vertex[plane]*ftimetick*fdriftvelocity); // in cm 
//     omega = asin(AC/sqrt(pow(AC,2)+pow(BC,2)));
//     omega = 180*omega/3.14; // in deg
//     //std::cout << " WireI1=" << wireI1 << " BI= " << BI << "    ThetaI = " << thetaI <<std::endl;
//        
//     if( (omega>(tan(slope[plane])-1.0*pi/180.))&&(omega<(tan(slope[plane])-1.0*pi/180.)) ){
//       fWire_last[plane] = wire;
//       fChannel_last[plane]=  channel;  // wire coordinate of vertex for each plane
//       fTime_last[plane] = time;
//     }
// 
// }   // end of HitIter loop
// 
// 
// 
// 
// 
//   
//   // Making the two lines in the two views
// for (unsigned int pl=0;pl<fNPlanes;pl++)
//   {
//   slope_wt[pl] = (fTime_last[pl]-fTime_vertex[pl])/(ftimetick*fdriftvelocity)/((fWire_last[pl]-fWire_vertex[pl])/ fMean_wire_pitch);
//   //slope_wt[0] = (fTime_last[0]-Time_C_wt)/(fWire_last[0]-Wire_vertexC_wt);
//  
//   intercept_wt[pl] = fTime_vertex[pl]/(ftimetick*fdriftvelocity) - fWire_vertex[pl]/(fMean_wire_pitch)*slope_wt[pl];
// //  intercept_wt[0] = Time_C_wt - Wire_vertexC_wt*slope_wt[0];
//   }
//   
//   return (void)0;
// }


int shwf::ShowerReco::Get3DaxisN(int iplane0,int iplane1){


 double l(0),m(0),n(0);
 double angle[3];
 // double Wire_vertex[3];
 //std::vector< double > angle;
  // Get Geometry
  art::ServiceHandle<geo::Geometry> geom;
 
unsigned int Cplane=0,Iplane=1;   // pretend collection and induction planes. "Collection" is the plane with the vertical angle equal to zero. If both are non zero collection is the one with the negative angle. 



//  slope[iplane0] = tan((fOmega_Mean[iplane0]*pi/180));
//  slope[iplane1] = tan((fOmega_Mean[iplane1]*pi/180));

//std::cout << " ---------- new calc, fOM, slope " << fOmega_Mean[iplane0] << " " << fOmega_Mean[iplane1] << " " <<  slope[iplane0] << " " << slope[iplane1] <<  std::endl;


 if(slope[iplane0]==0 || slope[iplane1]==0){
//    Phi   = 0;
//    Theta = 0;
   fPhiN_ang.push_back(-900);
  fThetaN_ang.push_back(-900);
    std::cout << " event parallell in one of the planes, exiting " << std::endl;
return 0;
  }
  


   // for now commented because this is already done in 3Daxis. This should be done in 2dvariables.

  //Wire_vertex[iplane0] = fWire_vertex[iplane0] *fMean_wire_pitch; // in cm 
  //Wire_vertex[iplane1] = fWire_vertex[iplane1] *fMean_wire_pitch; // in cm
  //Time_vertex[iplane0]  = fTime_vertex[iplane0]  *ftimetick*fdriftvelocity; // in cm

/////!!!! check this better.

 // intercept[iplane0] = fTime_vertex[iplane0]*ftimetick*fdriftvelocity - (slope[iplane0]*Wire_vertex[iplane0]);
 // intercept[iplane1] = fTime_vertex[iplane0]*ftimetick*fdriftvelocity - (slope[iplane1]*Wire_vertex[iplane1]);





//////////insert check for existence of planes.

  angle[iplane0]=geom->Plane(iplane0).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to 
  angle[iplane1]=geom->Plane(iplane1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to

  std::cout <<  " actual angles " << angle[iplane0]*180/pi << " " << angle[iplane1]*180/pi << std::endl;
  
  
 if(angle[iplane0]==0)   // first plane is at 0 degrees
       {
       Cplane=iplane0;
       Iplane=iplane1;
std::cout << "+++++ new calc first test case 1 angle[0]==0 "  << std::endl;
	}
else if(angle[iplane1]==0)  // second plane is at 0 degrees
       {
	Cplane=iplane1;
       Iplane=iplane0;
std::cout << "+++++ new calc first test case 2 angle[1]==0 "  << std::endl;
	}
else if(angle[iplane0]!=0 && angle[iplane1]!=0)  //both planes are at non zero degree - find the one with deg<0
      {
	if(angle[iplane1]<angle[iplane0])
          {Cplane=iplane1;
          Iplane=iplane0;
std::cout << "+++++ new calc first test case 3 angle[1]< angle[0] "  << std::endl;
          }
	else if(angle[iplane1]>angle[iplane0])
          {Cplane=iplane0;
          Iplane=iplane1;
std::cout << "+++++ new calc first test case 4 angle[1]> angle[0] "  << std::endl;
          }
        else
	  {
	//throw error - same plane.
           return -1;
	  }	

	}




//0 -1 factor depending on if one of the planes is vertical.
bool nfact=!(angle[Cplane]);


l = 1;

std::cout << "+++++ new calc first test c(I) ,I,C,I/C " << cos(angle[Iplane]) << " " <<  cos(angle[Cplane]) << " " << cos(angle[Iplane])/cos(angle[Cplane]) <<" nfact  " << nfact << " sines, I,C " << sin(angle[Iplane]) << " " << sin(angle[Cplane]) << std::endl;



  m = (1/(2*sin(angle[Iplane])))*((cos(angle[Iplane])/(slope[Cplane]*cos(angle[Cplane])))-(1/slope[Iplane]) +nfact*(  cos(angle[Iplane])/slope[Cplane]-1/slope[Iplane]  )     );


 /*double mm = (1/(2*sin(angle[Iplane])))*((1/slope[Cplane])-(1/slope[Iplane]));
double mmm = (1/(sin(angle[Iplane])))*((cos(angle[Iplane])/slope[Cplane])-(1/slope[Iplane]));
double mmmm = -(1/(sin(angle[Iplane])))*(((-1)*cos(angle[Iplane])/slope[Cplane])+(1/slope[Iplane]));
double mxm= (1/(sin(angle[Iplane])))*((cos(angle[Iplane])/(slope[Cplane]*cos(angle[Cplane])))-(1/slope[Iplane]));
double mcm=(1/(sin(angle[Iplane])))*( nfact*(  cos(angle[Iplane])/slope[Cplane]-1/slope[Iplane]  )     ); */
// std::cout << "+++++ new calc first test M vs. two non zero: " << m << " " << mm << std::endl;
// std::cout << "+++++ new calc first test M vs. one non zero (I=2): " << m << " " << mmm << " " << mmmm<< " " << mxm << " " <<mcm << std::endl;


  n = (1/(2*cos(angle[Cplane])))*((1/slope[Cplane])+(1/slope[Iplane])+nfact*((1/slope[Cplane])-(1/slope[Iplane])));
// double nn = (1/(2*cos(angle[Cplane])))*((1/slope[Cplane])+(1/slope[Iplane])); 
// 
// std::cout << "+++++ new calc first test N: " << n << " " << nn << std::endl;
// std::cout << "+++++ new calc first test N: " << n << " " << 1/slope[Cplane] << std::endl;
// 
//   std::cout << "________________m= " << m << std::endl;


  // Director angles
  fPhiN.push_back( atan(n/l));
  fThetaN.push_back( acos(m/(sqrt(pow(l,2)+pow(m,2)+pow(n,2)))) );


std::cout << "+++++ new calc first test angles tests, Phi, Thet: " << fPhiN[fPhiN.size()-1]*180/pi << " " << fThetaN[fThetaN.size()-1]*180/pi << std::endl;

  float Phi = fPhiN[fPhiN.size()-1]>0. ? (TMath::Pi()/2)-fPhiN[fPhiN.size()-1] : fabs(fPhiN[fPhiN.size()-1])-(TMath::Pi()/2) ; // solve the ambiguities due to tangent periodicity
  float Theta=0;
  if(Phi<0)Theta = (TMath::Pi()/2)-fThetaN[fThetaN.size()-1];
  if(Phi>0)Theta = fThetaN[fThetaN.size()-1]-(TMath::Pi()/2);


//std::cout << "+++++ new calc first test angles tests, after tangent corr, Phi, Thet: " << fPhiN[fPhiN.size()-1]*180/pi << " " << fThetaN[fThetaN.size()-1] << std::endl << std::endl;
 
  std::cout << " NPhi=" <<Phi*180/pi << "   Ntheta=" << Theta*180/pi <<std::endl;
  
  
    fh_phi_act->Fill(Phi*180/pi);

    fh_thet_act->Fill(Theta*180/pi);
  

  
  GetPitchLength(); //Get pitch of (two) wire planes (only for Argoneut)
   

  fPhiN_ang.push_back(Phi*180/pi);
  fThetaN_ang.push_back(Theta*180/pi);





  //slope[0] = tan((thetaC*pi/180));




return 0;

}







// int shwf::ShowerReco::Get3Daxis_coords(){
// 
//  art::ServiceHandle<geo::Geometry> geom;
// 
// int nvertices=0;
// 
// if(fNPlanes==2)
//    nvertices=1;
// if(fNPlanes==3)
//    nvertices=3;
// 
// fYvertex.resize(0);fZvertex.resize(0);
// fXvertex.resize(0);fZlast.resize(0);
// fYlast.resize(0);fXlast.resize(0);
// 
// //for(int iplane=0;nvertices<;iplane++)
// //   {
// double y,z;
// bool wires_cross = geom->ChannelsIntersect(fChannel_vertex[0],fChannel_vertex[1],y,z);
// fYvertex.push_back(y);fZvertex.push_back(z);
// wires_cross = geom->ChannelsIntersect(fChannel_last[0],fChannel_last[1],y,z);
// fYlast.push_back(y);fZlast.push_back(z);
// 
// fXvertex.push_back(fTime_vertex[0]*ftimetick*fdriftvelocity);
// fXlast.push_back(fTime_last[0]*ftimetick*fdriftvelocity ); 
// 
// TVector3 XYZ0;  // track origin or interaction vertex
// XYZ0.SetXYZ(fXlast[0]-fXvertex[0],fYlast[0]-fYvertex[0],fZlast[0]-fZvertex[0]);
// 
// std::cout << "^^^^^ channels " << fChannel_vertex[0] << " " << fChannel_vertex[1] << " " << fChannel_last[0] << " " << fChannel_last[1] << std::endl;
// std::cout << "^^^^^ X,Y,Z start " << fXvertex[0] << " " << fYvertex[0] << " " << fZvertex[0] << std::endl;
// std::cout << "^^^^^ X,Y,Z end " << fXlast[0] << " " << fYlast[0] << " " << fZlast[0] << std::endl;
// std::cout << "^^^^^ X,Y,Z diff " << fXlast[0]-fXvertex[0] << " " << fYlast[0]-fYvertex[0] << " " << fZlast[0]-fZvertex[0] << std::endl;
// std::cout << "^^^^^ angle calculation  theta: " << XYZ0.Theta() << " " << XYZ0.Phi() << std::endl; 
// 
// 
// 
// wires_cross = geom->ChannelsIntersect(fChannel_vertex[0],fChannel_vertex[2],y,z);
// fYvertex.push_back(y);fZvertex.push_back(z);
// wires_cross = geom->ChannelsIntersect(fChannel_last[0],fChannel_last[2],y,z);
// fYlast.push_back(y);fZlast.push_back(z);
// 
// fXvertex.push_back(fTime_vertex[0]*ftimetick*fdriftvelocity);
// fXlast.push_back(fTime_last[0]*ftimetick*fdriftvelocity ); 
// 
// TVector3 XYZ1;  // track origin or interaction vertex
// XYZ0.SetXYZ(fXlast[1]-fXvertex[1],fYlast[1]-fYvertex[1],fZlast[1]-fZvertex[1]);
// 
// 
// std::cout << "^^^^^ channels " << fChannel_vertex[0] << " " << fChannel_vertex[2] << " " << fChannel_last[0] << " " << fChannel_last[2] << std::endl;
// std::cout << "^^^^^ X,Y,Z start " << fXvertex[1] << " " << fYvertex[1] << " " << fZvertex[1] << std::endl;
// std::cout << "^^^^^ X,Y,Z end " << fXlast[1] << " " << fYlast[1] << " " << fZlast[1] << std::endl;
// std::cout << "^^^^^ X,Y,Z diff " << fXlast[1]-fXvertex[1] << " " << fYlast[1]-fYvertex[1] << " " << fZlast[1]-fZvertex[1] << std::endl;
// std::cout << "^^^^^ angle calculation  theta: " << XYZ1.Theta() << " " << XYZ1.Phi() << std::endl; 
// 
// 
// 
// wires_cross = geom->ChannelsIntersect(fChannel_vertex[1],fChannel_vertex[2],y,z);
// fYvertex.push_back(y);fZvertex.push_back(z);
// wires_cross = geom->ChannelsIntersect(fChannel_last[1],fChannel_last[2],y,z);
// fYlast.push_back(y);fZlast.push_back(z);
// 
// fXvertex.push_back(fTime_vertex[0]*ftimetick*fdriftvelocity);
// fXlast.push_back(fTime_last[0]*ftimetick*fdriftvelocity ); 
// 
// 
// TVector3 XYZ2;  // track origin or interaction vertex
// XYZ0.SetXYZ(fXlast[2]-fXvertex[2],fYlast[2]-fYvertex[2],fZlast[2]-fZvertex[2]);
// 
// std::cout << "^^^^^ channels " << fChannel_vertex[1] << " " << fChannel_vertex[2] << " " << fChannel_last[1] << " " << fChannel_last[2] << std::endl;
// std::cout << "^^^^^ X,Y,Z start " << fXvertex[2] << " " << fYvertex[2] << " " << fZvertex[2] << std::endl;
// std::cout << "^^^^^ X,Y,Z end " << fXlast[2] << " " << fYlast[2] << " " << fZlast[2] << std::endl;
// std::cout << "^^^^^ X,Y,Z diff " << fXlast[2]-fXvertex[2] << " " << fYlast[2]-fYvertex[2] << " " << fZlast[2]-fZvertex[2] << std::endl;
// std::cout << "^^^^^ angle calculation  theta: " << XYZ2.Theta() << " " << XYZ2.Phi() << std::endl; 
// 
// 
//  //  }
// 
// return 0;
// }







 
// int shwf::ShowerReco::Get3Daxis(float thetaI, float thetaC, float Wire_vertexI, float Wire_vertexC, float Time_vertex){
//   
//   //Get theta and phi (polar angles "direction of the shower")
// 
//   //float ftimetick = 0.198;    //time sample in microsec
//  
//   slope[1] = tan((thetaI*pi/180));
//   slope[0] = tan((thetaC*pi/180));
// 
// 
// std::cout << "$$$$$$$$$$$ slopes " <<  slope[0] << " " << slope[1] << "  " << thetaI <<" "<<thetaC << std::endl << std::endl;
// 
// /*  Wire_vertexI = Wire_vertexI *fMean_wire_pitch; // in cm 
//   Wire_vertexC = Wire_vertexC *fMean_wire_pitch; // in cm
//   Time_vertex  = Time_vertex  *ftimetick*fdriftvelocity; // in cm*/
// 
//   intercept[1] = Time_vertex - (slope[1]*Wire_vertexI);
//   intercept[0] = Time_vertex - (slope[0]*Wire_vertexC);
// 
//  
//   float l(0),m(0),n(0);
//   // Get Geometry
//   art::ServiceHandle<geo::Geometry> geom;
//  
//   float Angle = geom->Plane(0).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
//   std::cout << "Angle= " <<  Angle<<std::endl;
//   
//   float angle_rad = 30* pi /180;
//   angle_rad = Angle;
//   
//   l = 1;
//   m = (1/(2*sin(angle_rad)))*((1/slope[0])-(1/slope[1]));
//   n = (1/(2*cos(angle_rad)))*((1/slope[0])+(1/slope[1]));
//  
//   std::cout << "________________m= " << m << std::endl;
// 
//   // Director angles
//   fPhi   = atan(n/l);
//   fTheta = acos(m/(sqrt(pow(l,2)+pow(m,2)+pow(n,2))));
// 
//   float Phi = fPhi>0. ? (TMath::Pi()/2)-fPhi : fabs(fPhi)-(TMath::Pi()/2) ; // solve the ambiguities due to tangent periodicity
//   float Theta=0;
//   if(Phi>0)Theta = fTheta-(TMath::Pi()/2);
//   if(Phi<0)Theta = (TMath::Pi()/2)-fTheta;
// 
//   if(slope[0]==0 || slope[1]==0){
//     Phi   = 0;
//     Theta = 0;
//     fTheta = 0;
//     fPhi = 0;
//   }
//   
//   std::cout << " Phi=" <<Phi*180/pi << "   theta=" << Theta*180/pi <<std::endl;
//  
//   
//   
//     fh_phi_act->Fill(Phi*180/pi);
// 
//     fh_thet_act->Fill(Theta*180/pi);
//   
// 
//   
//   //GetPitchLength(); //Get pitch of (two) wire planes (only for Argoneut)
//    
// 
//   fPhi_ang=Phi*180/pi;
//   fTheta_ang=Theta*180/pi;
// 
//   
//   return 0;
// }


void shwf::ShowerReco::LongTransEnergy(art::PtrVector < recob::Hit> hitlist)
{
  // alogorithm for energy vs dx of the shower (roto-translation) COLLECTION VIEW
 // double  wire_cm, time_cm;
 // int loop_nrg = 0;

  
  double totCnrg = 0,totCnrg_corr =0 ; // tot enegry of the shower in collection
  //double CdEdx4cm = 0; // tot enegry of the shower in collection
  //int CdedxCounter = 0;
  art::ServiceHandle<geo::Geometry> geom;

  int channel;

//double CdEdx4cm=0,CdEdx6cm=0,CdEdx8cm=0,CdEdx10cm=0;

 
  double time;
  unsigned int wire=0,plane=0;

        
  for(unsigned int i=0;i<fNPlanes;i++)
  {
    	fTotChargeADC[i]=0;
	fTotChargeMeV[i]=0;
      //  fChannel_vertex.resize(fNPlanes,-1);
      //  fChannel_last.resize(fNPlanes,-1);

fChargeADC_4cm[i]=0;   //Initial charge in ADC/cm for each plane angle calculation 4cm
  fChargeMeV_4cm[i]=0;  //initial charge in MeV/cm for each angle calculation first 4cm
  fChargeADC_6cm[i]=0;   //Initial charge in ADC/cm for each angle calculation first 6cm
  fChargeMeV_6cm[i]=0;  //initial charge in MeV/cm for each angle calculation first 6cm
  fChargeADC_8cm[i]=0;   //Initial charge in ADC/cm for each angle calculationfirst 8cm
  fChargeMeV_8cm[i]=0;  //initial charge in MeV/cm for each angle calculation first 8cm
  fChargeADC_10cm[i]=0;   //Initial charge in ADC/cm for each angle calculation first 10cm
  fChargeMeV_10cm[i]=0;  //initial charge in MeV/cm for each angle calculation first 10cm

    
  }
  


  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
    //    if(time_C<1020)continue;
  //  wire_cm = wire * fMean_wire_pitch; //in cm
  //  time_cm = time *ftimetick*fdriftvelocity; //in cm
   
     // moving to polar coordinates
    
 
    double lifetime = exp(time*ftimetick/735.);
///////////// !!! To Be corrected significantly for plane difference and different values.

    totCnrg +=(theHit->Charge()*lifetime)/fNPitch[plane]; // Sum the energy of all the hits
    double dQdx_e = (theHit->Charge(true)*lifetime)/fNPitch[plane]/(calFactor*eCharge);  // in e/cm
///////////// !!! To be corrected for calFactor per plane? 
      totCnrg_corr += BirksCorrection(dQdx_e); 
   

double hlimit=30.4;  //arbitrary value to get first ~30 cm of the shower.    



double wdist=((wire-fWire_vertex[plane])*fNPitch[plane]);

if( (wdist<hlimit)&&(wdist>0.4)){ 
      //if(CdedxCounter<=10)

      
      
        double ch_adc=(theHit->Charge()*lifetime)/fNPitch[plane];
        double dQdx_e = (theHit->Charge(true)*lifetime)/fNPitch[plane]/(calFactor*eCharge);
        double Bcorr =  BirksCorrection(dQdx_e);
	
	//std::cout << " charge calc, plane " << plane << " " << fNPitch[plane] << " " << wdist << " " << ch_adc << " " << dQdx_e << " " << fChargeADC_4cm[plane] << " " << fChargeMeV_4cm[plane] << std::endl;
	
        // fill out for 4cm preshower
	  if(wdist<4.4)
             {fChargeADC_4cm[plane]+=ch_adc; 
              fChargeMeV_4cm[plane]+= Bcorr ; }
        // fill out for 6cm preshower
       	  if(wdist<6.4)
             {fChargeADC_6cm[plane]+=ch_adc;
             fChargeMeV_6cm[plane]+= Bcorr; }
        
         // fill out for 8cm preshower
       	  if(wdist<8.4)
             {fChargeADC_8cm[plane]+=ch_adc;
              fChargeMeV_8cm[plane]+= Bcorr; }
        
         // fill out for 10cm preshower
       	  if(wdist<10.4)
             {fChargeADC_10cm[plane]+=ch_adc;
              fChargeMeV_10cm[plane]+= Bcorr; }
         
        fDistribChargeADC[plane].push_back(ch_adc);  //vector with the first De/Dx points
	fDistribChargeMeV[plane].push_back(Bcorr);  //vector with the first De/Dx points
        fDistribChargeposition[plane].push_back(wdist);  //vector with the first De/Dx points' positions 


      //CdedxCounter++;
      //std::cout << " plane nr " << plane << " - hit-Vertex=" << ((wire_cm-fWire_vertex[plane]*fMean_wire_pitch)/fMean_wire_pitch)*fNPitch[plane] << " dedx|(MIPs/cm)=" << theHit->Charge()/fNPitch[plane] << std::endl;
      fh_dedx[plane]->Fill( wdist, ((((theHit->Charge()/lifetime)/fNPitch[plane])*10/7)/7.6)*6250*23.6*pow(10,-6) )  ;


    }
    
  }
  
 std::cout << " ------------ plane nr " << plane << " - hit-Vertex=" << std::endl;


       fTotChargeADC[plane]=totCnrg; 
       fTotChargeMeV[plane]=totCnrg_corr;  
  
  //  std::cout << "COLECTION - dEdx|4cm= " << CdEdx4cm<< "   TotCnrg= " << totCnrg << "  in GeV="  <<(((((totCnrg/7.3)*6300)*100)/60)*23.6)*pow(10,-9)<< std::endl;

}

//------------------------------------------------------------------------------------//  




// ******************************* //

// // Angular distribution of the energy of the shower - Collection view
// void shwf::ShowerReco::AngularDistribution(art::PtrVector < recob::Hit>  hitlist){
//   std::cout << "------ in angular distribution, n of hits " << hitlist.size() << std::endl;
//   int    loop = 0; // flag
//   art::ServiceHandle<geo::Geometry> geom;
//   double time;
//   unsigned int wire;
//   double BC,AC;
//   double omega;
//   unsigned int channel,plane;
// 
//   // this should changed on the loop on the cluster of the shower
//   for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
//     art::Ptr<recob::Hit> theHit = (*hitIter);
//     time = theHit->PeakTime();  
//     //time_C -= (presamplings+10.1);
//     art::Ptr<recob::Wire> theWire = theHit->Wire();
//     channel = theWire->RawDigit()->Channel();
//     geom->ChannelToWire(channel, tpc, plane, wire);
//     
//     //Here we should take GetVertex function to retrieve the vertex of the shower
//   //  if(loop==0){
//     //  fWire_vertex[plane] = wire;
//     //   fChannel_vertex[plane]=channel;
//     //   fTime_vertex[plane] = time;
//      // std::cout << "+++++++ setting vertex @ plane: "<< plane << " " << fWire_vertex[plane] << " " << fTime_vertex[plane] << std::endl;
//    // }
// 
//     //    std::cout << "VertexWireC= " << fWire_vertex[0] << "   VerTimeC= " << fTime_vertex[0] << std::endl;
//     
//     // moving to polar coordinate
//     BC = (wire - fWire_vertex[plane])*fMean_wire_pitch + fMean_wire_pitch; // in cm
//     AC = (time - fTime_vertex[plane])*ftimetick*0.158; //in cm 
//     omega = asin(  AC/sqrt(pow(AC,2)+pow(BC,2)) );
//     omega = 180*omega/3.14;
//     fh_theta[plane]->Fill(omega, theHit->Charge()); // Filling the histo (angle, energy of the hit)
//     fh_omega_evt[plane]->Fill(omega, theHit->Charge());
//     fh_omega_evt_reb[plane]->Fill(omega, theHit->Charge());
//     loop++; // flag counter
//   }
//   std::cout << "VertexWireC= " << fWire_vertex[plane] << "   VerTimeC= " << fTime_vertex[plane] << std::endl;
// 
//   return (void)0;
// }

// automatically called by Get3Daxis
void shwf::ShowerReco::GetPitchLength(){
  
  //Get pitch of 2 wire planes 
  //for generalization to n planes, different formulas for the geometrical transofrmation are needed (this is ArgoNeuT specific)
  // Get Geometry
 

  // TPC parameters
//   TString tpcName = geom->GetLArTPCVolumeName();
//   float Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
//   std::cout << "theta_input=" <<theta << " Phi=" << phi <<std::endl;  
//   
//   float angle_rad = (180+30)* pi /180;
//   fPitch[0] = fabs(fmean_wire_pitch/((cos(angle_rad)*sin(theta)*sin(phi))+(sin(angle_rad)*cos(theta))));
//   fPitch[1] = fabs(fmean_wire_pitch/((cos(angle_rad)*sin(theta)*sin(phi))-(sin(angle_rad)*cos(theta))));
// 
//   std::cout << "IND-Pitch[1]=" <<fPitch[1] << " COL-Pitch[0]=" << fPitch[0] <<std::endl;   
// 

//   
//   std::cout << "INDUCTION  - dEdx|4cm= " << IdEdx4cm<< "   " <<IdedxCounter<<std::endl;
//   std::cout << "COLLECTION - dEdx|4cm= " << CdEdx4cm<< "   " <<CdedxCounter<<std::endl;

//// Obtain pitch using the ProjectedLength method obtained from prong; 
 art::ServiceHandle<geo::Geometry> geom;

for(unsigned int pl=0;pl<fNPlanes;pl++)  
   {fNPitch[pl]=ProjectedLength(pl);
    std::cout << "++++++ calculating  -Pitch[" <<pl<<"] =" <<fNPitch[pl] << std::endl;

  } 	
//std::cout << std::endl; 

// fNPitch[1]=ProjectedLength(geom->Plane(1).View());

  return (void)0;
}


// /*****************************************************/
// void shwf::ShowerReco::GetVertex(art::Event& evt){
// 
//   // here we should GetThe verteces list
//   // then understand which one belongs to a shower
//   // getting the cluster that contain that vertex and pass everything to the main routine for the pghysical variables determination
//  // Read in the hough line List object(s).
// //     art::Handle< std::vector<recob::Cluster> > houghListHandle;
// //     evt.getByLabel(fHoughLineModuleLabel,houghListHandle);
// 
// art::Handle< std::vector<recob::Cluster> > lineListHandle;
//     evt.getByLabel(fVertexCLusterModuleLabel,lineListHandle);
// 
// 
//     
// 
//   /* // Read in the cluster List object(s).
//     art::Handle< std::vector<recob::Cluster> > clusterListHandle;
//     evt.getByLabel(fClusterModuleLabel,clusterListHandle);
//     // Read in the vertex Strength List object(s).
//     art::Handle< std::vector<recob::EndPoint2D> > vertexStrengthListHandle;
//     evt.getByLabel(fVertexStrengthModuleLabel,vertexStrengthListHandle);
//   */  
// 
// art::ServiceHandle<geo::Geometry> geo;
// 
// double wire_low[3]={10000,10000,10000},time_low[3]={100000,10000,10000};
// 
// // for(unsigned int ii = 0; ii < houghListHandle->size(); ++ii)
// for(unsigned int ii = 0; ii <lineListHandle->size(); ++ii)
//     {
// 
//     unsigned int local_wire_low;
//     double local_time_low;
// 
//       art::Ptr<recob::Cluster> hl(lineListHandle, ii);
//       
//       // Figure out which View the cluster belongs to 
//       
//       ///      int clPlane = cl->View()-1;
// 
//       art::PtrVector<recob::Hit> hitlist;
//       hitlist = hl->Hits();
//       hitlist.sort(shwf::SortByWire());
//       unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire
// 
//       art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin();
// 	c=(*a)->Wire()->RawDigit()->Channel(); 
// 	geo->ChannelToWire(c,t,p,local_wire_low);
//        local_time_low = (*a)->PeakTime();  
// 
// 
// 	std::cout << "@@@@@@@  planeview in hit list   " <<  geo->Plane(p,t).View() << " w_low "<<local_wire_low << " t_low " <<  local_time_low << std::endl; 
//         
//         if(local_wire_low < wire_low[geo->Plane(p,t).View()-1] && hitlist.size()>7 )
//              {
//              wire_low[geo->Plane(p,t).View()-1]=local_wire_low;
//              time_low[geo->Plane(p,t).View()-1]=local_time_low;
//               }
// 
// //          hitlist_all[geo->Plane(p,t).View()-1].push_back(*a);
// 
// 
// 
//     } // End loop on clusters.
// 
// 
// 
// for(int i=0;i<fNPlanes;i++)
// {std::cout << "@@@@@@@  ending values   " <<  i << " w_low "<<wire_low[i]*fMean_wire_pitch << " t_low " <<  time_low[i]*ftimetick*fdriftvelocity << std::endl;
// std::cout << "@@@@@@@  ending values   " <<  i << " w_low "<<wire_low[i] << " t_low " <<  time_low[i] << std::endl;
// fWire_vertex[i]=wire_low[i];
// fTime_vertex[i]=time_low[i];
// 
// }
//  
//   return (void)0;
// }



/*

void shwf::ShowerReco::GetVertexN(art::Event& evt){

 art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fMCGeneratorLabel,mctruthListHandle);




 art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);	
      mclist.push_back(mctparticle);
    } 




//for(int j = 0; j < mc->NParticles(); ++j){
  //    simb::MCParticle part(mc->GetParticle(j));


        //event_has_pi_plus=1;
      

//double vertex[3]={0,0,0};
//for( unsigned int i = 0; i < mclist.size(); ++i )
  //  art::Ptr<simb::MCTruth> mc(mclist[0]);

std::cout << "%%%%%%% mc size size,  "<<mclist.size() <<    std::endl;


    art::Ptr<simb::MCTruth> mc(mclist[0]);
    simb::MCParticle neut(mc->GetParticle(0));

if((neut.PdgCode()==22)&& neut.StatusCode()==1) //photon - need first electron.
     { 
     art::Handle< std::vector<sim::Particle> > parHandle;
     evt.getByLabel(fLarGeantlabel, parHandle);

     art::PtrVector<simb::MCParticle> pvec;
    int fpart=0;
    for(unsigned int i = 0; i < parHandle->size(); ++i){
      art::Ptr<simb::MCParticle> p(parHandle, i);      
      pvec.push_back(p);
      if(p->PdgCode() ==11 || p->PdgCode()==-11)
	  {
		fpart=i;
		break;
	  }	

    }


     std::cout << "%%%&&&&&&&&&& is PDG: " << pvec[fpart]->PdgCode() << " " << pvec[fpart]->TrackId() << std::endl;
      

    int trackid ;

	trackid =  mc->GetParticle(0).Daughter(0);
        std::cout << "####### NDaughters: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() <<  " "<<mc->GetParticle(1).NumberDaughters() <<  std::endl;

	for(int xx=0;xx<mc->GetParticle(0).NumberDaughters();xx++)
		{
      trackid =  mc->GetParticle(0).Daughter(xx);
        std::cout << "####### is PDG, trackid: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() << std::endl; 
              } 
	unsigned int jj;
      for(jj = 0; jj < pvec.size(); jj++) // Don't look below i.
	    {
	      if (trackid==pvec[jj]->TrackId())
		{
		   std::cout << "daughter particle "<<jj << " " << pvec[jj]->PdgCode() << std::endl; // get the pointer, 
			break; 
            
                 }            


    	   }
     neut=*pvec[fpart];

     } //end foton clause
//if((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1)
//NumberDaughters	(		 ) 	 const [inline]
//int simb::MCParticle::Daughter	(	const int 	i	 ) 	 const
 
 //  std::cout << "%%%%%%% particle size, partslist: "<< partslist << " "  <<  mc->NParticles() << std::endl; 
  int npart=0;
   //  while(&& npart < mc->NParticles() )
     //     {
                 std::cout << "%%%%%%%####### is PDG: "<< npart <<" " << neut.PdgCode() << std::endl; 
 	//	neut=mc->GetParticle(npart++);

       //   }       

 std::cout << "%%%%%%%####### after loop is PDG: "<< npart <<" " << neut.PdgCode() << std::endl; 
    //if((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1){
  
    
    xyz_vertex[0] =neut.Vx();
    xyz_vertex[1] =neut.Vy();
    xyz_vertex[2] =neut.Vz();
	
    std::cout<<"neut.Vx()= "<<neut.Vx()<<" ,y= "<<neut.Vy()<<" ,z= "<<neut.Vz()<<std::endl;
//if(((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1))
  //    break;

	
    
   art::ServiceHandle<geo::Geometry> geom;
art::ServiceHandle<util::LArProperties> larp;
double drifttick=(xyz_vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);

const double origin[3] = {0.};
for(unsigned int p=0;p<fNPlanes;p++)
{
double pos[3];
unsigned int  wirevertex, t;
geom->Plane(p).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
std::cout << "plane X positionp " << p << " " << pos[0] << std::endl;

pos[1]=xyz_vertex[1];
pos[2]=xyz_vertex[2];
 unsigned int channel2 = geom->NearestChannel(pos);
       geom->ChannelToWire(channel2,t,p,wirevertex); 
       
fWire_vertex[p]=wirevertex;
fTime_vertex[p]=drifttick+(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);
std::cout<<"wirevertex= "<<wirevertex<< " timevertex " << fTime_vertex[p] <<std::endl;
}



  return (void)0;
}*/





 void   shwf::ShowerReco::GetVertexAndAnglesFromCluster(art::Ptr< recob::Cluster > clust,unsigned int plane) // Get shower vertex and slopes.
{




slope[plane]=clust->dTdW()*(ftimetick*fdriftvelocity)/fMean_wire_pitch;  //convert to cm/cm units needed in the calculation
fWire_vertex[plane]=clust->StartPos()[0];
fTime_vertex[plane]=clust->StartPos()[1];

////////// insert detector offset

std::cout << "======= setting slope for view: " << plane << " " << slope[plane] << " " << fWire_vertex[plane] << " " << fTime_vertex[plane] << " " <<  fWire_vertex[plane]+50<< " "<< fTime_vertex[plane] + slope[plane]*(fWire_vertex[plane]+50) << std::endl;

}







/*****************************************************/
	
double shwf::ShowerReco::BirksCorrection(double dQdx_e){
  /// Correction for charge quenching using parameterization from 
  /// S.Amoruso et al., NIM A 523 (2004) 275
  // copied from CaloArgoItaliano.cxx
  double dQdx = dQdx_e;
  double dEdx;
  float A3t = 0.800; 
  float K3t = 0.0486; // in KV/cm*(g/cm^2)/MeV
  float rho = 1.388; // LAr density in g/cm^3
  double Wion = 23.6e-6 ; //   23.6 eV = 1e, Wion in MeV/e
  double Efield = 0.485;  // Electric Field in the drift region in KV/cm
  K3t = K3t/rho; // KV/MeV
  dEdx = dQdx/(A3t/Wion-K3t/Efield*dQdx); //MeV/cm 
  return dEdx;
  }



/***************************************************/
// taken and modified from recob::Prong           //


double shwf::ShowerReco::ProjectedLength(unsigned int plane ) 	 const
{
   art::ServiceHandle<geo::Geometry> geo;
    double pitch = -1.;
    if(geo->Plane(plane).View() == geo::kUnknown || geo->Plane(plane).View() == geo::k3D){
      mf::LogWarning("RecoBaseProng")<< "Warning Prong::ProjectedLength :  no Pitch foreseen for view "<<geo->Plane(plane).View();
      return pitch;
    }
    else{
 // temporary override of global angles to use the new calculations.
      double fTheta=fThetaN[0];
      double fPhi=fPhiN[0];  

     
      double pi=TMath::Pi();
      for(unsigned int i = 0; i < geo->Nplanes(); ++i){
        if(i == plane){
          double wirePitch = geo->WirePitch(0,1,i);
          double angleToVert =0.5*TMath::Pi() - geo->Plane(i).Wire(0).ThetaZ(false) ;

	  std::cout <<" %%%%%%%%%%  " << i << " angle " << angleToVert*180/pi << " " <<geo->Plane(i).Wire(0).ThetaZ(false)*180/pi <<" wirePitch " << wirePitch  <<std::endl;
           std::cout <<" %%%%%%%%%%  " << fTheta << " " << fPhi << std::endl;
	   
	   //(sin(angleToVert),cos(angleToVert)) is the direction perpendicular to wire
          //(fDCosStart[1],fDDosStart[2]) is the direction of prong in the y-z plane
          double cosgamma = TMath::Abs(TMath::Sin(angleToVert)*TMath::Cos(fTheta)+TMath::Cos(angleToVert)*TMath::Sin(fTheta)*TMath::Sin(fPhi));
          if (cosgamma>0) pitch = wirePitch/cosgamma;     
        } // end if the correct view
      } // end loop over planes
    } // end if a reasonable view

    return pitch;
  }


