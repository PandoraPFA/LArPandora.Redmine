////////////////////////////////////////////////////////////////////////
//
// \file ShowerAngleCluster.cxx
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
#include "Simulation/sim.h"
#include "ClusterFinder/ShowerAngleCluster.h"
#include "ClusterFinder/EndPointService.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"


#include "SimulationBase/simbase.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/summary.h"


// ***************** //


cluster::ShowerAngleCluster::ShowerAngleCluster(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
}


void cluster::ShowerAngleCluster::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel = pset.get< std::string >("ClusterModuleLabel");
   //fHoughLineModuleLabel=pset.get<std::string > ("HoughLineModuleLabel");
  fVertexCLusterModuleLabel=pset.get<std::string > ("VertexClusterModuleLabel");
  fMCGeneratorLabel=pset.get<std::string > ("MCGeneratorLabel");
  fLarGeantlabel=pset.get<std::string > ("LarGeantlabel");     
  fUseMCVertex=pset.get<int > ("UseMCVertex");
 

}

// ***************** //
cluster::ShowerAngleCluster::~ShowerAngleCluster()
{
}

namespace cluster {
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
void cluster::ShowerAngleCluster::beginJob()
{



  /** Get Geometry*/
  art::ServiceHandle<geo::Geometry> geo;
  fNPlanes = geo->Nplanes();
  fMean_wire_pitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  unsigned int tpc=0;
  
  art::ServiceHandle<util::DetectorProperties> detp;
  ftimetick=detp->SamplingRate()/1000.;
  
  std::cout << "------------ timetick? " << ftimetick << std::endl << std::endl << std::endl;

  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  /** Create Histos names*/
 
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  
  
 
  
  
  
  
  for(unsigned int i=0;i<fNPlanes;++i){


    //    sprintf(&tit_dedx[0],"fh_dedx_%.4i_%.4i_%i",i);
   
   int nwires=geo->Plane(i,tpc).Nwires();
   int ntimes=geo->DetHalfWidth(tpc)*2/(ftimetick*0.158);
//ftimetick      =  0.198; // time sample in us fdriftvelocity =  0.157;	
std::cout << "{{{ --- }}} init for plane" << i << " " << nwires << " " << ntimes << std::endl; 


    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_theta_%i",i);
    fh_theta[i] = tfs->make<TH1F>(tit_h_theta,"Theta distribution",720,-180., 180.);

    sprintf(&tit_h_theta[0],"charge distrib_%i",i);
    tgx[i]=tfs->make<TH2F>(tit_h_theta,"charge distribution per wires",nwires/8.,0, nwires*fMean_wire_pitch,ntimes/8.,0,ntimes*ftimetick*0.158);
   sprintf(&tit_h_theta[0],"hit distrib_%i",i);
    tgx2[i]=tfs->make<TH2F>(tit_h_theta,"Hit distribution per wires",nwires/8.,0, nwires*fMean_wire_pitch,ntimes/8.,0,ntimes*ftimetick*0.158);  

    std::cout << " ------ test hist widths " << ntimes*ftimetick*0.158 << " det half " << ntimes << std::endl; 

    linefit[i]=tfs->make<TF1>(Form("linefit_%d",i),"pol1",0,4000);
    linefit2[i]=tfs->make<TF1>(Form("linefit_2_%d",i),"pol1",0,4000);	
    
    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_omega_evt_%i",i);
    fh_omega_evt.push_back( tfs->make<TH1F>(tit_h_theta,"Theta distribution per event",720,-180., 180.) );

   
   sprintf(&tit_h_theta[0],"fh_omega_evt_reb_%i",i);
    fh_omega_evt_reb.push_back( tfs->make<TH1F>(tit_h_theta,"Theta distribution per event, rebinned",180,-180., 180.) );
   

   //fh_omega2_evt.push_back( new TH1F(tit_h_theta,"Theta distribution per event",720,-180., 180.) );    


    /**Histos for the angular distribution theta wire of the shower*/
    sprintf(&tit_h_theta[0],"ftheta_wire_%i",i);
    fh_theta_wt[i] = tfs->make<TH1F>(tit_h_theta_wt,"Theta wire distribution",720,-180., 180.);
 	}
  
  ftree_cluster =tfs->make<TTree>("ShowerAngleCluster","Results");/**All-knowing tree with reconstruction information*/
  
  
   ftree_cluster->Branch("run",&fRun,"run/I");
    ftree_cluster->Branch("subrun",&fSubRun,"subrun/I");
   ftree_cluster->Branch("event",&fEvent,"event/I");
   ftree_cluster->Branch("nplanes",&fNPlanes,"nplanes/I");
  
     ftree_cluster->Branch("mcpdg",&mcpdg,"mcpdg/I");
    ftree_cluster->Branch("mcenergy",&mcenergy,"mcenergy/D");
   
    ftree_cluster->Branch("mcphi",&mcphi,"mcphi/D");
    ftree_cluster->Branch("mctheta",&mctheta,"mctheta/D");
   
    
    
   ftree_cluster->Branch("wire_vertex","std::vector<unsigned int>", &fWire_vertex);
   ftree_cluster->Branch("time_vertex","std::vector<double>", &fTime_vertex);
   ftree_cluster->Branch("wire_last","std::vector<unsigned int>", &fWire_last);
   ftree_cluster->Branch("time_last","std::vector<double>", &fTime_last);

   ftree_cluster->Branch("test_wire_start","std::vector<double>", &test_wire_start);
   ftree_cluster->Branch("test_time_start","std::vector<double>", &test_time_start);
	



//   ftree_cluster->Branch("fitw_vertex","std::vector<double>", &wire_start);
//    ftree_cluster->Branch("fitt_vertex","std::vector<double>", &time_start);
   ftree_cluster->Branch("fitw_last","std::vector<double>", &wire_end);
   ftree_cluster->Branch("fitt_last","std::vector<double>", &time_end); 

 
  
  
    ftree_cluster->Branch("xyz_vertex","std::vector<double>", &xyz_vertex);
    ftree_cluster->Branch("xyz_vertex_fit","std::vector<double>", &xyz_vertex_fit);

 // hioni->Branch("Nhit",  Nhit ,"Nhit[Npart]/I");
// this should be temporary - until the omega is sorted out.
  // ftree_cluster->Branch("fh_omega2_evt","std::vector<TH1F *>", &fh_omega2_evt);



    ftree_cluster->Branch("omega_2d","std::vector<double>", &fOmega_Mean);
    ftree_cluster->Branch("omega_2d_RMS","std::vector<double>", &fOmega_RMS);

    ftree_cluster->Branch("omega_2d_line","std::vector<double>", &fOmega_Mean_line);
    ftree_cluster->Branch("omega_2d_RMS_line","std::vector<double>", &fOmega_RMS_line);

    ftree_cluster->Branch("omega_2d_reb","std::vector<double>", &fOmega_Mean_reb);
    ftree_cluster->Branch("omega_2d_reb_RMS","std::vector<double>", &fOmega_RMS_reb);
    ftree_cluster->Branch("omega_2d_mean","std::vector<double>", &fOmega_Mean_Mean);

ftree_cluster->Branch("slope","std::vector<double>", &slope);
ftree_cluster->Branch("lineslope","std::vector<double>", &lineslope);
ftree_cluster->Branch("calcslope","std::vector<double>", &calcslope);

ftree_cluster->Branch("RMS_wire","std::vector<double>", &fRMS_wire);
ftree_cluster->Branch("RMS_time","std::vector<double>", &fRMS_time);

ftree_cluster->Branch("Chisq","std::vector<double>", &fChisq);
ftree_cluster->Branch("minwir","std::vector<double>", &fminwir);
ftree_cluster->Branch("maxwir","std::vector<double>", &fmaxwir);
ftree_cluster->Branch("mintime","std::vector<double>", &fmintime);

ftree_cluster->Branch("maxtime","std::vector<double>", &fmaxtime);
ftree_cluster->Branch("correlation","std::vector<double>", &fcorrelation);
ftree_cluster->Branch("covariance","std::vector<double>", &fcovariance);



   ftree_cluster->Branch("Eventangleposition","std::vector<std::vector<double>>",&fSingleEvtAngle);
ftree_cluster->Branch("Eventanglepositionval","std::vector<std::vector<double>>",&fSingleEvtAngleVal);

  // ftree_cluster->Branch("fslope_2d"," std::vector<double>", &fSlope_2d);
  // ftree_cluster->Branch("fintercept_2d","std::vector<double>", &fIntercept_2d);
//   
ftree_cluster->Branch("ShowerPosition2D","std::vector<std::vector<double>>",&fShowerPosition2D);
ftree_cluster->Branch("ShowerWidthProfile2D","std::vector<std::vector<double>>",&fShowerWidthProfile2D);
ftree_cluster->Branch("ShowerChargeProfile2D","std::vector<std::vector<double>>",&fShowerChargeProfile2D);
  


}

// ***************** //
void cluster::ShowerAngleCluster::produce(art::Event& evt)
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
  
  std::cout << "------------ timetick? " << ftimetick << std::endl << std::endl << std::endl;
  
//     	fWire_vertex.resize(0);  // wire coordinate of vertex for each plane
//     	fTime_vertex.resize(0);  // time coordinate of vertex for each plane
// 	fWire_last.resize(0);  // wire coordinate of vertex for each plane
//     	fTime_last.resize(0);  // time coordinate of vertex for each plane
// 
//         fOmega_Mean.resize(0);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
//         fOmega_RMS.resize(0);;     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm
// 	
//         fOmega_Mean_reb.resize(0);    // Mean value of the 2D angular Rebinned by 4
//         fOmega_RMS_reb.resize(0);     // RMS of the 2D angular distribution  Rebinned by 4
//         fOmega_Mean_Mean.resize(0);    // Mean value of the 2D angular use mean instead of maximum
//         
// 
//         fOmega_wt_Mean.resize(0);; // Mean value of the angular distribution (1=Ind - 0=Coll) wire,time
//         fOmega_wt_RMS.resize(0);;  // RMS of the angular distribution  (1=Ind - 0=Coll) wire,time
//         fChannel_vertex.resize(0);  // wire coordinate of vertex for each plane
//          fChannel_last.resize(0);  // wire coordinate of vertex for each plane



	//fPitch.resize(0);  // Pitch calculated the old way
    	

 
fSingleEvtAngle.resize(fNPlanes); 
fSingleEvtAngleVal.resize(fNPlanes); 

 fShowerWidthProfile2D.resize(fNPlanes); ;  // vector to show the plane shower Width distribution 
 fShowerChargeProfile2D.resize(fNPlanes); ;  //vector to show the plane shower Charge distribution
 fShowerPosition2D.resize(fNPlanes); ;  //vector to store the positions of hit values stored in the previous two vectors.


 for(unsigned int ii=0;ii<fNPlanes;ii++)
  {   
  fSingleEvtAngle[ii].resize(180); 
  fSingleEvtAngleVal[ii].resize(180); 
  }


      // fPitch.resize(fNPlanes); 
    	 
	fWire_vertex.resize(fNPlanes);
	fTime_vertex.resize(fNPlanes);
        fWire_last.resize(fNPlanes);
	fTime_last.resize(fNPlanes);
 	fChannel_vertex.resize(fNPlanes);
        fChannel_last.resize(fNPlanes);

	xyz_vertex.resize(fNPlanes);
	xyz_vertex_fit.resize(fNPlanes);	

       test_wire_start.resize(fNPlanes);
       test_time_start.resize(fNPlanes);

   wire_start.resize(fNPlanes);wire_end.resize(fNPlanes);;
   time_start.resize(fNPlanes);time_end.resize(fNPlanes);;

	slope.resize(fNPlanes);
	lineslope.resize(fNPlanes);
	calcslope.resize(fNPlanes);

        fOmega_Mean.resize(fNPlanes);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
        fOmega_RMS.resize(fNPlanes);     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm

        fOmega_Mean_line.resize(fNPlanes);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
        fOmega_RMS_line.resize(fNPlanes);     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm

        fOmega_wt_Mean.resize(fNPlanes); // Mean value of the angular distribution (1=Ind - 0=Coll) wire,time
        fOmega_wt_RMS.resize(fNPlanes);  // RMS of the angular distribution  (1=Ind - 0=Coll) wire,time
        fOmega_Mean_reb.resize(fNPlanes);    // Mean value of the 2D angular Rebinned by 4
        fOmega_RMS_reb.resize(fNPlanes);     // RMS of the 2D angular distribution  Rebinned by 4
        fOmega_Mean_Mean.resize(fNPlanes);    // Mean value of the 2D angular use mean instead of maximum
 fRMS_wire.resize(fNPlanes);
 fRMS_time.resize(fNPlanes);
 fChisq.resize(fNPlanes);
 fminwir.resize(fNPlanes);
 fmaxwir.resize(fNPlanes);
 fmintime.resize(fNPlanes);
 fmaxtime.resize(fNPlanes);
 fcorrelation.resize(fNPlanes);
 fcovariance.resize(fNPlanes);         




  /**Get Clusters*/
  std::cout << "************ What I'm getting out " << fClusterModuleLabel << " " << std::endl;
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);


   std::vector< art::PtrVector < recob::Hit> > hitlist_all;
   //art::PtrVector < recob::Hit> hitlistInd;
   hitlist_all.resize(fNPlanes);

 
    //std::auto_ptr<std::vector<recob::Shower> > Shower3DVector(new std::vector<recob::Shower>);

   
   art::ServiceHandle<cluster::EndPointService>  endpoints;
    art::PtrVector<recob::Cluster> clusters;


  unsigned int clmaxsize[3]={0.};
  int clmax[3]={0.}; 
   //first loop on clusters to find the biggest one in each planes
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii){

      art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      
  
     
      
      art::PtrVector<recob::Hit> hitlist;
      hitlist = cl->Hits();
      unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire
      art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin();
      c=(*a)->Wire()->RawDigit()->Channel(); 
      geo->ChannelToWire(c,t,p,w);
      
      if(hitlist.size()>clmaxsize[p]){
	clmaxsize[p]=hitlist.size();
	clmax[p]=ii;
      }
      
      
    } // end loop determining largest clusters.
      
      
      for(unsigned int ii=0;ii<fNPlanes;ii++)
      {
	std::cout << "cluster size " << clmaxsize[ii] << " nclust " << clmax[ii] << std::endl;
       art::Ptr<recob::Cluster> cl(clusterListHandle, clmax[ii]);
       clusters.push_back(cl);
      }
      
      
      for(unsigned int ii = 0; ii < clusters.size(); ++ii){

      
     
      
      art::PtrVector<recob::Hit> hitlist;
      hitlist = clusters[ii]->Hits();
      
      unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire

      art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin();
      c=(*a)->Wire()->RawDigit()->Channel(); 
      geo->ChannelToWire(c,t,p,w);
      wire_start[p]=w;

//       a = hitlist.end();
//       c=(*a)->Wire()->RawDigit()->Channel(); 
//       geo->ChannelToWire(c,t,p,w);
//       wire_end[p]=w;	

      for(art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin(); a != hitlist.end();  a++) //loop over cluster hits
      {
	c=(*a)->Wire()->RawDigit()->Channel(); 
	geo->ChannelToWire(c,t,p,w);

        //geo->Plane(i).View()

	

          hitlist_all[p].push_back(*a);

// 	if(geo->Plane(p,t).SignalType() == geo::kCollection)
// 	  {	  
// 	    hitlistCol.push_back(*a);
// 
// 	  }
// 	else if (geo->Plane(p,t).SignalType() == geo::kInduction)
// 	  {
// 	    hitlistInd.push_back(*a);
// 	  } 
      }

    //  hitlistInd.sort(cluster::SortByWire());
     // hitlistCol.sort(cluster::SortByWire());


    } // End loop on clusters.


  std::cout << " ---- endpoint coordinates " << std::endl;

// for(unsigned int ii=0;ii<fNPlanes;ii++)
//   {
//   
//   }
//       std::vector< recob::EndPoint2D > vtxcol;
//        endpoints->EndPoint(clusters,vtxcol);
//        for(unsigned int yy=0;yy<vtxcol.size();yy++)
// 	  {
// 	   std::cout << "       ++++ plane, elem " << vtxcol[yy].View() << " " << yy << " pos " << vtxcol[yy].WireNum() << " " << vtxcol[yy].DriftTime() << " strength " << vtxcol[yy].Strength() << std::endl;
// 	    
// 	  }




 fRun = evt.id().run();
 fSubRun = evt.id().subRun();
 fEvent = evt.id().event();



   // GetVertex(evt);
if(fUseMCVertex)
    GetVertexN(evt);


for(unsigned int i=0;i<fNPlanes;i++)
   AngularDistribution(hitlist_all[i]); // 2D Direction of the shower in consecutive planes

 Find2DStartPoints( hitlist_all);
  

//move to 2Dstartpoints
for(unsigned int p=0;p<fNPlanes;p++)
{fWire_vertex[p]=wire_start[p];
fTime_vertex[p]=time_start[p];
}

   for(unsigned int i=0;i<fNPlanes;i++)
      {
       hitlist_all[i].sort(cluster::SortByWire());
      fh_omega_evt[i]->Reset();
      fh_omega_evt_reb[i]->Reset();
     
       FitAngularDistributions(hitlist_all[i]);  
      Get2DVariables(hitlist_all[i]);
	}

 
  //AngularDistribution(hitlistCol); // 2D Direction of the shower Collection
              // Fit of 2d distributions (for both planes)
  
  std::cout << "######## in main loop " << fOmega_Mean[0] << " " <<  fOmega_Mean[1] << std::endl;
  
art::PtrVector<recob::Hit> clusterHits;



// make an art::PtrVector of the clusters
std::auto_ptr<std::vector<recob::Cluster> > ShowerAngleCluster(new std::vector<recob::Cluster>);
unsigned int i=0;

for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
{
 clusterHits.clear();
	int maxlength=0,maxpos=0;
 	for( i = 0; i < clusterListHandle->size(); ++i){
   		art::Ptr<recob::Cluster> prod(clusterListHandle, i);

		

      art::PtrVector<recob::Hit> hitlist;
      hitlist = prod->Hits();
      
      unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire

      art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin();
      c=(*a)->Wire()->RawDigit()->Channel(); 
      geo->ChannelToWire(c,t,p,w);                

// does not work due to plane reordering.
//if((prod->View()-1)!=iplane)
//			continue;

	if((p)!=iplane)
			continue;

	//	std::cout << "----- looping on clusters " << i << " " << prod->View()-1 << " " <<prod->Hits().size() << std::endl;
	//	std::cout << "----- lenght and pos " << i << " " << maxlength << " "  <<maxpos << std::endl;

   		if(  (prod->Hits().size() > (unsigned int)maxlength) )
		{
        		maxpos=i;
        		maxlength=prod->Hits().size();
          //        std::cout << "++++ condition ok " << std::endl << std::endl;
		}
 	}

art::Ptr<recob::Cluster> prod(clusterListHandle, maxpos);

for(unsigned int ii=0;ii<prod->Hits().size();ii++)
  clusterHits.push_back(prod->Hits()[ii]);


int wirevert=prod->StartPos()[0];
double timevert=prod->StartPos()[1];

int wireend=wire_end[iplane];
double timeend=time_end[iplane];


if(fUseMCVertex)
{
wirevert=fWire_vertex[iplane];
timevert=fTime_vertex[iplane];
}

recob::Cluster temp(clusterHits,wirevert, prod->SigmaStartPos()[0],timevert, prod->SigmaStartPos()[1], wireend, prod->SigmaEndPos()[0],timeend, prod->SigmaEndPos()[1], slope[iplane], slope[iplane]*0.05, lineslope[iplane],lineinterc[iplane], iplane);

std::cout << "######## in plane loop filling clusters " << std::endl; 

ShowerAngleCluster->push_back(temp);
}



  /**Fill the output tree with all information */
    ftree_cluster->Fill();


  evt.put(ShowerAngleCluster);

}


// ******************************* //

// Angular distribution of the energy of the shower - Collection view
void cluster::ShowerAngleCluster::AngularDistribution(art::PtrVector < recob::Hit>  hitlist){
  std::cout << "------ in angular distribution, n of hits " << hitlist.size() << std::endl;
  int    loop = 0; // flag
  art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire;
  unsigned int channel,plane;


 art::Ptr<recob::Hit> theHit = (*hitlist.begin());
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
  
    
unsigned int minwire=wire,maxwire=0;;
double mintime=99999,maxtime=0.;

	tgx[plane]->Reset();
        tgx2[plane]->Reset();
   	//tgx[plane]->Set(hitlist.size());

  // this should changed on the loop on the cluster of the shower
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
    
   maxwire=wire;   

   if(time>maxtime)
	maxtime=time;

if(time<mintime)
	mintime=time;

 

  

    loop++; // flag counter
  }
  std::cout << "VertexWireC= " << fWire_vertex[plane] << "   VerTimeC= " << fTime_vertex[plane] << std::endl;
 


  int nbinsx= (maxwire+20-minwire+20)*fMean_wire_pitch;  // nbins to have 
  int nbinsy= (maxtime+200-mintime+200)*ftimetick*0.158;  // nbins to have 

 
  tgx[plane]->SetBins(nbinsx,(minwire-20)*fMean_wire_pitch,(maxwire+20)*fMean_wire_pitch,nbinsy,(mintime-200)*ftimetick*0.158,(maxtime+200)*ftimetick*0.158);
  tgx2[plane]->SetBins(nbinsx,(minwire-20)*fMean_wire_pitch,(maxwire+20)*fMean_wire_pitch,nbinsy,(mintime-200)*ftimetick*0.158,(maxtime+200)*ftimetick*0.158);
 
  
  
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
  
   tgx[plane]->Fill((double)wire*fMean_wire_pitch,time*ftimetick*0.158,theHit->Charge());
   tgx2[plane]->Fill((double)wire*fMean_wire_pitch,time*ftimetick*0.158);		

  }
  
double w_bar=0,t_bar=0;
int nhits=0;

for(int i=0;i<tgx2[plane]->GetNbinsX();i++)
	for(int j=0;j<tgx2[plane]->GetNbinsY();j++)
		{
		if(tgx2[plane]->GetBinContent(i,j)<=3)
			tgx2[plane]->SetBinContent(i,j,0);
		else
		      {
			w_bar+=tgx2[plane]->GetXaxis()->GetBinCenter(i);
			t_bar+=tgx2[plane]->GetYaxis()->GetBinCenter(j);
			nhits++;
			}	
		}

 

w_bar/=nhits;
t_bar/=nhits;

double Sxx=0,Syy=0,Sxy=0;

for(int i=0;i<tgx2[plane]->GetNbinsX();i++)
	for(int j=0;j<tgx2[plane]->GetNbinsY();j++)
		{
		if(tgx2[plane]->GetBinContent(i,j)>3)
			{Sxx+=(tgx2[plane]->GetXaxis()->GetBinCenter(i)-w_bar)*(tgx2[plane]->GetXaxis()->GetBinCenter(i)-w_bar);
			Syy+=(tgx2[plane]->GetYaxis()->GetBinCenter(j)-t_bar)*(tgx2[plane]->GetYaxis()->GetBinCenter(j)-t_bar);
			Sxy+=(tgx2[plane]->GetYaxis()->GetBinCenter(j)-t_bar)*(tgx2[plane]->GetXaxis()->GetBinCenter(i)-w_bar);;
			}
		}


 


 
  tgx[plane]->Fit(Form("linefit_%d",plane),"QMRNCFrob=0.8");
  tgx2[plane]->Fit(Form("linefit_2_%d",plane),"QMRNCFrob=0.95");


std::cout << "{{{-----}}}  histo stats: rms w,t " << tgx[plane]->GetRMS(1) << " " << tgx[plane]->GetRMS(2) << " chisq " << linefit[plane]->GetChisquare()/linefit[plane]->GetNDF() << " max, min wires and times " <<
minwire << " " <<maxwire << " " <<  mintime << " " << maxtime << std::endl;


fRMS_wire[plane]=tgx[plane]->GetRMS(1);
fRMS_time[plane]=tgx[plane]->GetRMS(2);
fChisq[plane]=linefit[plane]->GetChisquare()/linefit[plane]->GetNDF();
fminwir[plane]=minwire;
fmaxwir[plane]=maxwire;
fmintime[plane]=mintime;
fmaxtime[plane]= maxtime;
fcorrelation[plane]=tgx[plane]->GetCorrelationFactor();
fcovariance[plane]=tgx[plane]->GetCovariance();

  return (void)0;
}








// ***************** //
void cluster::ShowerAngleCluster::FitAngularDistributions(art::PtrVector < recob::Hit> hitlist){
  /** Fit function of the angular distribution (cm,cm)*/
 // art::ServiceHandle<geo::Geometry> geo;
  //unsigned int planes = geo->Nplanes();
  //TF1 *gau = new TF1("gaus","gaus",-60, 60);
 std::cout << "------ in angular distribution, n of hits " << hitlist.size() << std::endl;

  art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire;
  double BC,AC;
  double omega;
  unsigned int channel,iplane,plane;

art::Ptr<recob::Hit> theHit = (*hitlist.begin());
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, iplane, wire);



   	//tgx[plane]->Set(hitlist.size());

  // this should changed on the loop on the cluster of the shower
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
  
      BC = ((double)wire - fWire_vertex[plane])*fMean_wire_pitch; // in cm
      AC = ((double)time - fTime_vertex[plane])*ftimetick*0.158; //in cm 
      omega = asin(  AC/sqrt(pow(AC,2)+pow(BC,2)) );
 
      if(BC<0)  // for the time being. Will check if it works for AC<0
	  { 
	  if(AC!=0)
	  omega= AC/fabs(AC)*pi-omega;  //subtract when negative, add when positive
	  else    
	  omega=pi;
	  } 

      omega = 180*omega/3.14;
      fh_theta[plane]->Fill(omega, theHit->Charge()); // Filling the histo (angle, energy of the hit)
      fh_omega_evt[plane]->Fill(omega, theHit->Charge());

  
      fh_omega_evt_reb[plane]->Fill(omega, theHit->Charge());

  }
  //for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane){
  

    fOmega_Mean[iplane] =
    fh_omega_evt[iplane]->GetBinCenter(fh_omega_evt[iplane]->GetMaximumBin());// Mean value of the fit
    fOmega_RMS[iplane] = fh_omega_evt[iplane]->GetRMS(); // RMS of the fit of the angular distribution in deg

    fOmega_Mean_reb[iplane]= fh_omega_evt_reb[iplane]->GetBinCenter(fh_omega_evt_reb[iplane]->GetMaximumBin());// Mean value of the fit
    fOmega_RMS[iplane] = fh_omega_evt_reb[iplane]->GetRMS(); // RMS of the fit of the angular distribution in deg
    fOmega_Mean_Mean[iplane]= fh_omega_evt[iplane]->GetMean();// Mean value of the;    // Mean value of the 2D angular use mean instead of maximum
    
std::cout << "########## intermediate angles, plane: " << iplane << " stand, _w reb, mean " << fOmega_Mean[iplane] << " " << fOmega_Mean_reb[iplane] << " " << fOmega_Mean_Mean[iplane] << std::endl;


fOmega_Mean_line[iplane]=atan(linefit[iplane]->GetParameter(1));


for(int i=0;i<180;i++)
{fSingleEvtAngleVal[iplane][i]=fh_omega_evt_reb[iplane]->GetBinContent(i);
//fSingleEvtAngle[iplane][i]=fh_omega_evt[iplane]->GetBinContent(i);
fSingleEvtAngle[iplane][i]=(double)i*2-180;
}
  //}
//  double  Low_th  = fOmega_Mean[iplane]-(alpha*fOmega_RMS[iplane]);
//  double  High_th = fOmega_Mean[iplane]+(alpha*fOmega_RMS[iplane]);


slope[iplane] = tan((fOmega_Mean[iplane]*pi/180))*fMean_wire_pitch/(ftimetick*fdriftvelocity);



calcslope[iplane]=linefit2[iplane]->GetParameter(1);

std::cout << " ((------stand slope and slope from hits only ----- )) " << slope[iplane] << " " << calcslope[iplane] << "  "<< std::endl;

}















void cluster::ShowerAngleCluster::GetVertexN(art::Event& evt){

 art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fMCGeneratorLabel,mctruthListHandle);




 art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);	
      mclist.push_back(mctparticle);
    } 


std::cout << "%%%%%%% mc size size,  "<<mclist.size() <<    std::endl;


    art::Ptr<simb::MCTruth> mc(mclist[0]);
    simb::MCParticle neut(mc->GetParticle(0));

    mcpdg=neut.PdgCode();
    mcenergy=neut.P();  
    
    if (neut.P()){
      double lep_dcosx_truth = neut.Px()/neut.P();
      double lep_dcosy_truth = neut.Py()/neut.P();
      double lep_dcosz_truth = neut.Pz()/neut.P();

     std::cout << "-----  cx,cy,cz " << lep_dcosx_truth << " " << lep_dcosy_truth << " " << lep_dcosz_truth << std::endl;


mcphi=  (lep_dcosx_truth == 0.0 && lep_dcosz_truth == 0.0) ? 0.0 : TMath::ATan2(lep_dcosx_truth,lep_dcosz_truth);
mctheta= (lep_dcosx_truth == 0.0 && lep_dcosy_truth == 0.0 && lep_dcosz_truth == 0.0) ? 0.0 : TMath::Pi()*0.5-TMath::ATan2(sqrt(lep_dcosx_truth*lep_dcosx_truth + lep_dcosz_truth*lep_dcosz_truth),lep_dcosy_truth);




mcphi=180*mcphi/TMath::Pi();
mctheta= 180*mctheta/TMath::Pi();
    std::cout << "-----  phi, theta " <<  mcphi << " " << mctheta << std::endl;

    }

    
    
    
    
// if((neut.PdgCode()==22)&& neut.StatusCode()==1) //photon - need first electron.
//      { 
//      art::Handle< std::vector<sim::Particle> > parHandle;
//      evt.getByLabel(fLarGeantlabel, parHandle);
// 
//      art::PtrVector<simb::MCParticle> pvec;
//     int fpart=0;
//     for(unsigned int i = 0; i < parHandle->size(); ++i){
//       art::Ptr<simb::MCParticle> p(parHandle, i);      
//       pvec.push_back(p);
//       if(p->PdgCode() ==11 || p->PdgCode()==-11)
// 	  {
// 		fpart=i;
// 		break;
// 	  }	
// 
//     }
// 
// 
//      std::cout << "%%%&&&&&&&&&& is PDG: " << pvec[fpart]->PdgCode() << " " << pvec[fpart]->TrackId() << std::endl;
//       
// 
//     int trackid ;
// 
// 	trackid =  mc->GetParticle(0).Daughter(0);
//         std::cout << "####### NDaughters: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() <<  " "<<mc->GetParticle(1).NumberDaughters() <<  std::endl;
// 
// 	for(int xx=0;xx<mc->GetParticle(0).NumberDaughters();xx++)
// 		{
//       trackid =  mc->GetParticle(0).Daughter(xx);
//         std::cout << "####### is PDG, trackid: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() << std::endl; 
//               } 
// 	unsigned int jj;
//       for(jj = 0; jj < pvec.size(); jj++) // Don't look below i.
// 	    {
// 	      if (trackid==pvec[jj]->TrackId())
// 		{
// 		   std::cout << "daughter particle "<<jj << " " << pvec[jj]->PdgCode() << std::endl; // get the pointer, 
// 			break; 
//             
//                  }            
// 
// 
//     	   }
//      neut=*pvec[fpart];
// 
//      } //end foton clause
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
   

double drifttick=(xyz_vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);

const double origin[3] = {0.};
for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
{
double pos[3];
unsigned int  wirevertex, t;
unsigned int p;
geom->Plane(iplane).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
std::cout << "plane X positionp " << iplane << " " << pos[0] << std::endl;

pos[1]=xyz_vertex[1];
pos[2]=xyz_vertex[2];
 unsigned int channel2 = geom->NearestChannel(pos);
       geom->ChannelToWire(channel2,t,p,wirevertex); 
       
if(iplane!=p)
	{std::cout << " error - planes don't match " << iplane << " " << p << std::endl;
	return ;
	}

fWire_vertex[p]=wirevertex;
fTime_vertex[p]=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);
std::cout<<"wirevertex= "<<wirevertex<< " timevertex " << fTime_vertex[p] << " correction "<< (pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick) << " " << pos[0] <<std::endl;

}



  return (void)0;
}




//int cluster::ShowerAngleCluster::Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt){
void cluster::ShowerAngleCluster::Get2DVariables(art::PtrVector < recob::Hit> hitlist) {  

 
  // only needed for drawing the axis of the shower in the event display
  
 
 // double omega_sh, wire_cm, time_cm;


// double AC, BC, omega; 


//double minlength={10000};

///////// enter hits? or just wire coordinates?! 




//ok. now have, presumable start and endpoint positions for cluster. 




  return (void)0;
}






void   cluster::ShowerAngleCluster::Find2DStartPoints(std::vector< art::PtrVector < recob::Hit> > hitlist_all)
{


//// find for which planes the correlation factor is largest (and preferably larger than 0.6)

std::vector< int > best_planes;

for (unsigned int iplane=0;iplane<fNPlanes;iplane++)
    {
     if(fabs(tgx[iplane]->GetCorrelationFactor()) > 0.6 )	
	best_planes.push_back(iplane);
     std::cout << " correlation factors in 0.6 search" << tgx[iplane]->GetCorrelationFactor() << std::endl;
    }

 
 for (unsigned int iplane=0;iplane<fNPlanes;iplane++)
    {
     if(fabs(tgx[iplane]->GetCorrelationFactor()) <= 0.6 && fabs(tgx[iplane]->GetCorrelationFactor()) > 0.5 && tgx[iplane]->GetRMS(2) > 0. && tgx[iplane]->GetRMS(1) /tgx[iplane]->GetRMS(2)>2. )	
	best_planes.push_back(iplane);
     //std::cout << " correlation factors in 0.6 search" << tgx[iplane]->GetCorrelationFactor() << std::endl;
    }
 
 
    
//// Find which plane has the highest correlation factor, and if there is only one, add another with larger correl. factor
unsigned int used_plane=999;
while(best_planes.size()<2)
	{
	double mincorr=0;
	int maxplane=0;
	for (unsigned int iplane=0;iplane<fNPlanes;iplane++)
		{
		   std::cout << " in search loop search" << iplane << " used " << used_plane << " max " << maxplane << " mincorr " << mincorr << " " << tgx[iplane]->GetCorrelationFactor() << " boole ipl!=used_plane " << (iplane!=used_plane) << std::endl;
		if(fabs(tgx[iplane]->GetCorrelationFactor()) <= 0.6 && fabs(tgx[iplane]->GetCorrelationFactor()) > mincorr && iplane!=used_plane )	
			{maxplane=iplane;
			mincorr=fabs(tgx[iplane]->GetCorrelationFactor());
			
			}
		}

		std::cout << "pushing back " << maxplane << std::endl;
		used_plane=maxplane;  // to cut out redundancy.
	best_planes.push_back(maxplane);

	}

/////test values:
for(unsigned int ii=0;ii<best_planes.size();ii++)
	std::cout << "======+++==== " << best_planes[ii] << " " << tgx[best_planes[ii]]->GetCorrelationFactor() << std::endl; 


for (unsigned int iplane=0;iplane<fNPlanes;iplane++)
{
 std::cout << "other parameters: " <<  std::endl;
    std::cout << "RMS:: " <<  tgx[iplane]->GetRMS(1) << " " << tgx[iplane]->GetRMS(2) << " " << tgx[iplane]->GetRMS(1) /tgx[iplane]->GetRMS(2) << std::endl;
    TH1D *PX=tgx[iplane]->ProjectionX();
    std::cout << " nbins " << PX->GetNbinsX() << std:: endl;
    int maxix=-1,minix=-1;
    for(int i=0;i<PX->GetNbinsX();i++){
     if(PX->GetBinContent(i)>0){
       minix=i;
       break;
     }
    }
    for(int i=PX->GetNbinsX();i>0;i--){
     if(PX->GetBinContent(i)>0){
       maxix=i;
       break;
     }
    }
     
        
    TH1D *PY=tgx[iplane]->ProjectionY();
       std::cout << " nbins " << PY->GetNbinsX() << std:: endl;
    int maxiy=-1,miniy=-1;
    for(int i=0;i<PY->GetNbinsX();i++){
     if(PY->GetBinContent(i)>0){
       miniy=i;
       break;
     }
    }
    for(int i=PY->GetNbinsX();i>0;i--){
     if(PY->GetBinContent(i)>0){
       maxiy=i;
       break;
     }
    }
    std::cout << "max,min x,y:: " << PX->GetBinCenter(minix) << " "<< PX->GetBinCenter(maxix) << " " << PX->GetBinCenter(maxix) - PX->GetBinCenter(minix) << " y: " << PY->GetBinCenter(miniy) << " "<< PY->GetBinCenter(maxiy) << " " << PY->GetBinCenter(maxiy) - PY->GetBinCenter(miniy) << " "<<(PX->GetBinCenter(maxix) - PX->GetBinCenter(minix))/(PY->GetBinCenter(maxiy) - PY->GetBinCenter(miniy)) << std::endl;
    
    std::cout << " chisq " << linefit[iplane]->GetChisquare()<< " " << linefit[iplane]->GetNDF() << " "<< linefit[iplane]->GetChisquare()/linefit[iplane]->GetNDF() << std::endl;
 
}

 
 std::cout << std::endl << std::endl;
//// find the hits near the beginning for the best planes.

 art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire,plane;

double a,c;
double wlst,wlend;
double tlst,tlend;

double wire_bar=0,time_bar=0;
int nhits=0;
 unsigned int channel;

// loop on selected planes
for (unsigned int ii=0;ii<best_planes.size();ii++)
{

unsigned int iplane=best_planes[ii];

//get first wire of plane (should be sorted by wire number)
art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist_all[iplane].begin();
channel = (*hitIter)->Wire()->RawDigit()->Channel();
geom->ChannelToWire(channel, tpc, plane, wire);


//error checking:
if(iplane!=plane)
	{
	std::cout << " error: planes mismatch  " << iplane << " "<< plane << std::endl;
	return;
	}

//get paramters of the straight line fit.
a=linefit[iplane]->GetParameter(1);
lineslope[iplane]=a*fMean_wire_pitch/(ftimetick*fdriftvelocity);
 c=linefit[iplane]->GetParameter(0);
lineinterc[iplane]=c/(ftimetick*fdriftvelocity);

//get slope of lines orthogonal to those found crossing the shower.
 double aprim=0;
      
	if(a)	
	{
	aprim=-1./a;
	}

  std::cout << "========= line params, plane: a,c " << plane << " " << a << " " << slope[iplane] << " " << c << std::endl;




// start loop to find the extreme intercepts for points in the cluster - i.e. presumably the the lines going therough the first and last hit of the shower.

double extreme_intercept_end=-999999;
double extreme_intercept_start=999999;

double extr_wire_pos=0,extr_time_pos=0;

int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.

if(a>0)
	{
	multiplier=1;
	}
else if(a<0)
	{
	multiplier=-1;
	}

for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist_all[iplane].begin(); hitIter != hitlist_all[iplane].end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
 
   wire_bar+=wire;
   time_bar+=time;	
   nhits++;
    
      if(a)	
	{
// if a> 0 looking for larger intercept. if a<0 looking for smaller intercept as endpoint.
       

	double intercept=time*ftimetick*0.158-aprim*(double)wire*fMean_wire_pitch;
 //	double wire_on_line=(intercept - c)/(a-aprim);
   //     double time_on_line=a*wire_on_line+c;  

      std::cout << "[[----]] intercepts check, wire,time " << wire<< " " <<time << " int " << intercept << " end:  " << multiplier*extreme_intercept_end << " start: " << multiplier*extreme_intercept_start << std::endl;

	if(multiplier*intercept > extreme_intercept_end ) 	
		{
		extreme_intercept_end=intercept;
		//wire_end[plane]=wire;
		}

	if(multiplier*intercept < extreme_intercept_start ) 	
		{
		extreme_intercept_start=intercept;
	        extr_wire_pos=wire;extr_time_pos=time;
		//wire_start[plane]=wire;
		}  

	}

}   // end of first HitIter loop, at this point we should have the extreme intercepts 


std::cout << ":::::::::: extreme wire and time,start " << extr_wire_pos << " " << extr_time_pos << std::endl; 

wire_bar/=nhits;
time_bar/=nhits;

wlst=(extreme_intercept_start - c)/(a-aprim);  // in cm
tlst=(a*wlst+c);
//wlst=(extreme_intercept_start - c)/(a-aprim)/fMean_wire_pitch; // in wire number

wlend=(extreme_intercept_end - c)/(a-aprim);   //in cm
tlend=(a*wlend+c);   // in cm
//wlend=(extreme_intercept_end - c)/(a-aprim)/fMean_wire_pitch;

double Sxx=0,Sxy=0,Syy=0;

std::cout << "^^^^^^^^^ a^prim + max and min intercept " << aprim << " " << extreme_intercept_end << " " << extreme_intercept_start << std::endl;

double min_length_from_start=99999,min_length_from_end=99999;
int wire_online_end=(extreme_intercept_end - c)/(a-aprim); // in cm
int wire_online_begin=(extreme_intercept_start - c)/(a-aprim); // in cm
double time_online_end=a*wire_online_end+c; // in cm
double time_online_begin=a*wire_online_begin+c;  // in cm


std::cout << " :::::::: wire_online begin point " << wire_online_begin << " " << time_online_begin << std::endl;

//third loop to find first and last points (in theory)
for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist_all[iplane].begin(); hitIter != hitlist_all[iplane].end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);

	//double intercept=time*-aprim*(double)wire;
	double intercept=time*ftimetick*0.158-aprim*(double)wire*fMean_wire_pitch;
 	double wire_on_line=(intercept - c)/(a-aprim); //in cm
        double time_on_line=a*wire_on_line+c;     //in cm

//	double dist=TMath::Sqrt( pow(wire_on_line-(double)wire,2)+pow(time_on_line-time,2) );
//	double dist_from_start=TMath::Sqrt( pow(wire_on_line-wlst,2)+pow(time_on_line-tlst,2) );


    //double dist_begin=TMath::Sqrt( pow(wire_online_begin-(int)wire,2)+pow(time_online_begin-time,2) );
    double dist_begin_mod=TMath::Sqrt( pow(((double)wire_online_begin-(double)wire*fMean_wire_pitch),2)+pow((time_online_begin-time*0.158*ftimetick),2) );	


    Sxx+=(wire-wire_bar)*(wire-wire_bar);
    Syy+=(time-time_bar)*(time-time_bar);
    Sxy+=(time-time_bar)*(wire-wire_bar);

    if(dist_begin_mod<min_length_from_start)
	{
	wire_start[plane]=wire;
	time_start[plane]=time;
	min_length_from_start=dist_begin_mod;
	}	

   double dist_end=TMath::Sqrt( pow(((double)wire_online_end-wire*fMean_wire_pitch),2)+pow((time_online_end-time*0.158*ftimetick),2) );

    if(dist_end<min_length_from_end)
	{
	wire_end[plane]=wire;
	min_length_from_end=dist_end;
	}


//std::cout << "((----- )) point at " << wire << " " << time << "dist_b,_e : " << dist_begin << " " << dist_begin_mod << " " <<   min_length_from_start  << " ot dist "<< dist << " " << dist_from_start <<std::endl;


	fShowerPosition2D[plane].push_back(TMath::Sqrt( pow(wire_on_line-wlst,2)+pow(time_on_line-tlst,2) ));  
	fShowerWidthProfile2D[plane].push_back(TMath::Sqrt( pow(wire_on_line-(double)wire,2)+pow(time_on_line-time,2) ));
	fShowerChargeProfile2D[plane].push_back(theHit->Charge()); 
	



} // end of second hit loop. At this point we should have the wire and time closest to the beginning of the shower.


//calculate the first and last cluster points on the through line:
 wlst=wire_start[plane]*fMean_wire_pitch; // in cm
 tlst=a*wlst+c;   // in cm
 //time_start[plane]=tlst; 

 wlend=wire_end[plane]*fMean_wire_pitch;   // temporary - will need to get last wire coordinate for each plane and hitlist.
 tlend=a*wlend+c;  // in cm

time_end[plane]=tlend/(ftimetick*0.158);

//std::cout << "======== start and end positions for plane" << plane << " " << wlst << " " << tlst << " " << wlend << " " << tlend << std::endl; 



//std::cout << ":::::::: Sxx, Syy, Sxy " << Sxx << " " << Syy << " " << Sxy << " slope " << Sxy/Sxx << " corr vars " << Sxx/nhits<< " " << Syy/nhits <<" " << Sxy/nhits<< std::endl;


} // end of  loop on selected planes


for(unsigned int ii=0;ii<best_planes.size();ii++)
	std::cout << " ----++++ determined start wire points per planes " << best_planes[ii] << " "  << wire_start[best_planes[ii]] << " " << time_start[best_planes[ii]] <<std::endl;  


//sort the best planes in order of correlation factor


///////////////// Only do the following part if there are 3 planes.


if(fNPlanes>=3)
{


	//first sort
	if(tgx[best_planes[0]]->GetCorrelationFactor()<tgx[best_planes[1]]->GetCorrelationFactor())
		std::swap(best_planes[0],best_planes[1]);


	//second sort
	if(best_planes.size()>2 && (tgx[best_planes[1]]->GetCorrelationFactor()<tgx[best_planes[2]]->GetCorrelationFactor()))
		std::swap(best_planes[1],best_planes[2]);


	// sanity check - see if times of the points found are close enough.
	
	const double origin[3] = {0.};
	std::vector <std::vector <  double > > position;

	for(unsigned int xx=0;xx<best_planes.size();xx++)
		{
		double pos1[3];
	
		geom->Plane(best_planes[xx]).LocalToWorld(origin, pos1);
		std::vector <double > pos2;
		pos2.push_back(pos1[0]);
		pos2.push_back(pos1[1]);
		pos2.push_back(pos1[2]);
		position.push_back(pos2);
		}

	
	//loop to check time discrepancies between points found
	for(unsigned int xx=0;xx<best_planes.size()-1;xx++)  // best_planes.size()-1 because we want to always find a next plane
		{
		for(unsigned int yy=xx+1;yy<best_planes.size();yy++)  
			{
		std::cout << "**** difference between planes X position, planes: "<< best_planes[xx] << " "<< best_planes[yy] << " "<< position[xx][0] << " " << position[yy][0] << " " << fabs(position[xx][0]-position[yy][0]) << " " << " " << time_start[best_planes[xx]] << " " << time_start[best_planes[yy]] << " " << fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) << std::endl;
			
			if(fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) > 1.5*fabs(position[xx][0]-position[yy][0]))
				{
				std::cout << " !!!! Something's wrong in the wire determination " << fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) << " " << 1.5*fabs(position[xx][0]-position[yy][0]) << std::endl;
				}	
			} // end inner yy loop
		} // end outer xx loop

	
	

	if((fabs(time_start[best_planes[0]]-time_start[best_planes[1]]) > 1.5*fabs(position[0][0]-position[1][0])) && best_planes.size() > 2) //time discrepancy of best correlation factor pair and we have a choice.
		{std::cout << " time discrepancy of best correlation factor pair 0 and 1 " << time_start[best_planes[0]] << " " << time_start[best_planes[1]] << " "<< position[0][0] << " " << position[1][0] << std::endl;
		if(!(fabs(time_start[best_planes[0]]-time_start[best_planes[2]]) > 2.5*fabs(position[0][0]-position[2][0]))) //0,1 is bad but 0,2 is ok:
			{
			std::cout << " but 0,2 is ok " << time_start[best_planes[0]] << " " << time_start[best_planes[2]] << " "<< position[0][0] << " " << position[2][0] << std::endl;
			std::swap(best_planes[1],best_planes[2]);
			}
		else   //0,1 is not ok and 0,2 is not ok
			{
			std::cout << " 0,1 and 0,2 is not ok " << std::endl;
			if(!(fabs(time_start[best_planes[1]]-time_start[best_planes[2]]) > 2.5*fabs(position[1][0]-position[2][0]))) //0,1 and 0,2 is bad but 1,2 is ok.
				{
				std::cout << " but 1,2 is ok " << time_start[best_planes[1]] << " " << time_start[best_planes[2]] << " "<< position[1][0] << " " << position[2][0] << std::endl;
				std::swap(best_planes[0],best_planes[1]);
				std::swap(best_planes[1],best_planes[2]);  // shift zero to last position.
				}
			}

		}


	for(unsigned int ii=0;ii<best_planes.size();ii++)
	std::cout << " ----++++ determined start wire points per planes " << best_planes[ii] << " "  << wire_start[best_planes[ii]] << " " << time_start[best_planes[ii]]/(ftimetick*0.158) <<std::endl;  


	// Assuming there is no problem ( and we found the best pair that comes close in time )
	// we try to get the Y and Z coordinates for the start of the shower. 
	int chan1=geom->PlaneWireToChannel(best_planes[0],wire_start[best_planes[0]], 0);
	int chan2=geom->PlaneWireToChannel(best_planes[1],wire_start[best_planes[1]], 0);

	double y,z;
	bool wires_cross = geom->ChannelsIntersect(chan1,chan2,y,z);

	
	xyz_vertex_fit[1]=y;
	xyz_vertex_fit[2]=z;
	xyz_vertex_fit[0]=time_start[best_planes[0]]*0.158*ftimetick+position[0][0];


	std::cout << ":::::: found x,y,z vertex " << wires_cross << " " << xyz_vertex_fit[0] << " " << y << " " << z << std::endl;


	// assume some condition - probably that correlation factor is too small. Then project the found vertex in x,y,z into wire, time coordinates in the last plane.

	double pos[3];
	unsigned int  wirevertex, t;
	unsigned int worst_plane=2;
	if(best_planes.size()>=3)
		worst_plane=best_planes[2];
	else  //find plane that has bad correlation factor. We know, that at least two are there.
		{
		std::cout << "bplane size <3 " << best_planes.size() << std::endl;
		for(unsigned int jj=0;jj<fNPlanes;jj++)
			{
			bool exist_flag=false;		
			for(unsigned int kk=0;kk<best_planes.size();kk++)			
				{if(jj==(unsigned int)best_planes[kk])
					exist_flag=true;
				std::cout << " jj,kk, true or false " << jj << " " << kk << std::endl;	
				}

			if(!exist_flag)  // adding non_existing flag
				{
				worst_plane=jj;
				std::cout << "setting worst plane to " << jj << std::endl;
				break;
				}
			}
		}	

	geom->Plane(worst_plane).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
	std::cout << "plane X positionp " << worst_plane << " " << pos[0] << std::endl;

	
	
	///////////////////////////////////////
	// geometry test:
// 	double width  = 2.*geom->TPC(0).HalfWidth();  //notice the geometry gives the 1/2 width, so multiply by 2
// 	double height = 2.*geom->TPC(0).HalfHeight(); //notice the geometry gives the 1/2 height, so multiply by 2
// 	double length =    geom->TPC(0).Length();     //notice the geometry gives the total length
// 	
// 	
// 	std::cout << "-------- height " << height << " " << length << " " <<width << std::endl;
// 	
// 	//display first and last wires:
// 	
// 	for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	    { std::cout << " +++ pl" << iplane << geom->Plane(iplane).Nwires()  << std::endl;
// 	      double wire1_Start[3]={0},wire1_End[3]={0};
// 	    geom->WireEndPoints(0,iplane,0,wire1_Start,wire1_End);
// 	      std::cout << "++++++++ wire positions: w nr " << 0 << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	
// 	      geom->WireEndPoints(0,iplane,geom->Plane(iplane).Nwires()-1,wire1_Start,wire1_End);
// 	      std::cout << "++++++++ wire positions: w nr " << geom->Plane(iplane).Nwires()-1 << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	
// 	      
// 	
// 	    }
	
	// display all wire positions
	
// 	for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	{
// 	  std::cout << " +++ pl" << iplane << geom->Plane(iplane).Nwires()  << std::endl;
// 	for(unsigned int wire=0;wire<geom->Plane(iplane).Nwires();wire+=10)
// 	{
// 	  double wire1_Start[3]={0},wire1_End[3]={0};
// 	 geom->WireEndPoints(0,iplane,wire,wire1_Start,wire1_End);
// 	//std::cout << "++++++++ wire positions: w nr " << wire << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	 
// 	}
// 	std::cout << std::endl << std::endl << std::endl;
// 	}
// 	
// 	
// 	
// 	for(double iy=-60;iy<60;iy+=15)
// 	{
// 	  // for(double iz=length-500;iz<length-100;iz+=10)
// 	  for(double iz=1000;iz<1100;iz+=50)
// 	  {
// 	   unsigned int wirev[3]; 
// 	     std::cout << " ---- plane for y,z: " << iy << " " << iz <<" ";
// 	  for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	  {unsigned int p;
// 	  geom->Plane(iplane).LocalToWorld(origin, pos);
// 	  pos[1]=iy;
// 	  pos[2]=iz;
// 	    unsigned int channel2 = geom->NearestChannel(pos);
// 	  geom->ChannelToWire(channel2,t,p,wirev[iplane]); 
// 	 std::cout << " p: " << iplane << " " << wirev[iplane];
// 	  } // end of wire finding for positions
// 	   std::cout << std::endl;  
// 	  
// 	   // second loop to test intersect:
// 	   
// 	    std::cout << " ---- y,z for planes : " << std::endl; 
// 	   for(unsigned int iplane=0;iplane<fNPlanes-1;iplane++)
// 	      {
// 	      for(unsigned int iplane2=iplane+1;iplane2<fNPlanes;iplane2++)
// 		{
// 		int chan1=geom->PlaneWireToChannel(iplane,wirev[iplane], 0);
// 		int chan2=geom->PlaneWireToChannel(iplane2,wirev[iplane2], 0);
// 
// 		double y,z;
// 		bool wires_cross = geom->ChannelsIntersect(chan1,chan2,y,z);
// 	        
// 	       std::cout << iplane << " " << iplane2 << " " << y << " " << z << std::endl;
// 	      
// 	      wires_cross = geom->ChannelsIntersect(chan2,chan1,y,z);
// 	        
// 	       std::cout << " inverse " << iplane2 << " " << iplane << " " << y << " " << z << std::endl;
// 	      
// 	      
// 		}
// 	      }
// 	  
// 	 }
// 	}
	 
	//// end geometry test 
	/////////////////////////////////////////////////////////////////////////////////////////////// 
	 
	pos[1]=xyz_vertex_fit[1];
	pos[2]=xyz_vertex_fit[2];
 	unsigned int channel2 = geom->NearestChannel(pos);
       	geom->ChannelToWire(channel2,t,worst_plane,wirevertex); 


	art::ServiceHandle<util::LArProperties> larp;
	art::ServiceHandle<util::DetectorProperties> detp;
	double drifttick=(xyz_vertex_fit[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);


	double timestart=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick)+detp->TriggerOffset();
	std::cout << " worst plane " << worst_plane <<" wirevertex= "<<wirevertex<< " timevertex " << timestart << " correction " << detp->TriggerOffset() << " " << (pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick) << " "<<pos[0] <<std::endl;


	double min_dist=999999.;

	for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist_all[worst_plane].begin(); hitIter != hitlist_all[worst_plane].end();  hitIter++){
    		art::Ptr<recob::Hit> theHit = (*hitIter);
    		time = theHit->PeakTime() ;  
    		//time_C -= (presamplings+10.1);
    		unsigned int plane;
    		art::Ptr<recob::Wire> theWire = theHit->Wire();
    		channel = theWire->RawDigit()->Channel();
    		geom->ChannelToWire(channel, tpc, plane, wire);

	
    		double dist_begin=TMath::Sqrt( pow((double)((int)wirevertex-(int)wire)*fMean_wire_pitch,2)+pow((timestart-time)*0.158*ftimetick,2) );	

		//std::cout << "=== min_dist " << wire << " " << time <<" " << dist_begin << " " << pow((double)((int)wirevertex-(int)wire)*fMean_wire_pitch,2) << " " << ((int)wirevertex-(int)wire)*fMean_wire_pitch << " " << min_dist << std::endl; 
		
		if(dist_begin<min_dist)
			{
			min_dist=dist_begin;
			wire_start[worst_plane]=wire;
			time_start[worst_plane]=time;

			}


	} // end loop on hits.


} // end big if(fNPlanes >= 3)



std::cout << " final wire, time vertices for all planes:  " << std::endl;  

for(unsigned int ii=0;ii<fNPlanes;ii++)
	std::cout << ii << " " << wire_start[ii] << " " << time_start[ii] << " " << geom->Plane(ii,0).SignalType() << " " << geo::kCollection << " " << geo::kInduction  << std::endl;
	




}







