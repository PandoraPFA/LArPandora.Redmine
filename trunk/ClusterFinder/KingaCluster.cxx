////////////////////////////////////////////////////////////////////////
//
// KingaCluster class
//
// kinga.partyka@yale.edu
//

////////////////////////////////////////////////////////////////////////

#include "KingaCluster.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// ROOT includes
#include <TCanvas.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <dirent.h>

#include "TTree.h"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

 
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

cluster::KingaCluster::KingaCluster(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel       (pset.get< std::string >("DBScanModuleLabel")),
  fEndPoint2DModuleLabel        (pset.get< std::string >("EndPoint2DModuleLabel"))
  
{
  produces< std::vector<recob::Cluster> >();
}

cluster::KingaCluster::~KingaCluster()
{

}

//-------------------------------------------------
void cluster::KingaCluster::beginJob(){
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<geo::Geometry> geo;
  //unsigned int planes = geo->Nplanes();
 
 
 
  fh_theta_ind= tfs->make<TH1F>("fh_theta_ind","theta angle in degrees, Induction Plane", 180,-180 ,180  );
  fh_theta_coll= tfs->make<TH1F>("fh_theta_coll","theta angle in degrees, Collection Plane", 180,-180 ,180  );
  fh_theta_ind_2D= tfs->make<TH1F>("fh_theta_ind_2D","theta angle in degrees, Induction Plane", 180,-180 ,180  );
  fh_theta_coll_2D= tfs->make<TH1F>("fh_theta_coll_2D","theta angle in degrees, Collection Plane", 180,-180 ,180  );
 fh_theta_ind_Area= tfs->make<TH1F>("fh_theta_ind_Area","Hit Area vs theta angle in degrees, Induction Plane", 180,-180 ,180  );
  fh_theta_coll_Area= tfs->make<TH1F>("fh_theta_coll_Area","Hit Area vs theta angle in degrees, Collection Plane", 180,-180 ,180  );

Hit_Area_Ind= tfs->make<TH1F>("Hit_Area_Ind","Hit Area, Induction Plane", 100,0 ,1  );
Hit_Area_Coll= tfs->make<TH1F>("Hit_Area_Coll","Hit Area, Collection Plane", 100,0 ,1  );

}

//-----------------------------------------------------------------

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

void cluster::KingaCluster::produce(art::Event& evt)
{
  

std::cout<<"In KingaCluster::produce(art::Event& evt)"<<std::endl;
std::cout<<" Working on ";
std::cout << "Run: " << evt.run();
std::cout << " Event: " << evt.id().event() << std::endl;
  fpeaks_found=1;
  fMC=0; 
 //int RunNo= evt.run();
 //int EventNo=evt.id().event();
 art::ServiceHandle<util::LArProperties> larp;
// double electronlifetime=larp->ElectronLifetime();
 
ftime_vertex.clear();
fwire_vertex.clear();
ftime_vertex_reco.clear();
fwire_vertex_reco.clear();



if (evt.isRealData()) 
    {
      std::cout<<" YOU ARE WORKING WITH DATA !!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    //let's see if this Run has been scanned and thus we can find vertex info for it:
    // DIR *pDIR=0;
//     struct dirent *entry=0;
//     const char path[60]="/argoneut/data/simplescan_data/simplescan_text/";
//   char path_[60];
//   strcpy(path_,path);
//     if(pDIR==opendir(path))
//     {
//        while(entry==readdir(pDIR))
//        {
//          if(strcmp(entry->d_name,".")!=0 && strcmp(entry->d_name, "..")!=0)
//          {
//           std::string file=entry->d_name;
//           std::cout<<"*******file is "<<file<<std::endl;
//          std::string::size_type pos_end=file.rfind(".");
//          std::string no_string;
//          no_string=file.substr(pos_end-3,3);
//          std::cout<<"no_string= "<<no_string<<std::endl;
//          std::istringstream stream(no_string);
//          int No;
//          stream>>No;
//            if(RunNo==No)
//            {
//             std::cout<<"RUN NO "<<RunNo<<" was scanned, we can find a vertex info :) "<<std::endl;
//             //Now, open the file and find the vertex info in it:
//             //.........................................
//             strcat(path_,entry->d_name);
//             std::string k;
//             int run,event,blah1,blah2,blah3,blah4,time_ind,time_coll,w_ind,w_coll;
//             std::ifstream ScannedFile;
//             ScannedFile.open(path_,std::ios::in);
//             if(!ScannedFile.is_open()){std::cout<<" Couldn't open file named: "<<entry->d_name<<std::endl;}
//             if(ScannedFile.is_open()){std::cout<<" openED file named: "<<entry->d_name<<std::endl;} 
//              
//              
//             while(getline(ScannedFile,k))
//             {
//             std::istringstream ins;
//             ins.clear();
//             ins.str(k);
//             ins>>run>>event>>blah1>>blah2>>blah3>>blah4>>time_ind>>time_coll>>w_ind>>w_coll;
//             //std::cout<<run<<" "<<event<<std::endl;
//             if(EventNo==event && RunNo==run)
//             {
//             ftime_vertex.push_back(time_ind);
//             ftime_vertex.push_back(time_coll);
//             fwire_vertex.push_back(w_ind);
//             fwire_vertex.push_back(w_coll);
//             
//             std::cout<<"GOT VERTEX INFO FROM THE SCANNED FILE:"<<std::endl;
//             std::cout<<"("<<time_ind<<" , "<<w_ind<<" )"<<std::endl;
//             std::cout<<"("<<time_coll<<" , "<<w_coll<<" )"<<std::endl;
//             
//             
//             break;
//             }//event=event
//             
//             
//             
//             
//             
//             }//while getline
//             
//             
//             
//             
//            //.........................................
//             ScannedFile.close();
//             break; 
//          
//            } //if the run was scanned
//            
//            
//            
//          }
//        
//        
//        } //while it loops through all the entries
//         closedir(pDIR);
//     
//     
//     }//if can open directory
//     else{ std::cout<<" THIS RUN HASN'T BEEN SCANNED! SORRY! "<<std::endl;}
//     
    } //if realData




else {

std::cout<<" YOU ARE WORKING WITH MC !!!!!!!!!!!!!!!!!!!!! "<<std::endl;
 fMC=1;
 fGenieGenModuleLabel="generator";

art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
    
    
    
    
for( unsigned int i = 0; i < mclist.size(); ++i ){

    art::Ptr<simb::MCTruth> mc(mclist[i]);

    simb::MCParticle neut(mc->GetParticle(i));

    // std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
    MCvertex[0] =neut.Vx();
    MCvertex[1] =neut.Vy();
    MCvertex[2] =neut.Vz();
    std::cout<<"MCvertex[0]= "<<MCvertex[0]<<std::endl;
    std::cout<<"driftvelocity= "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())<<std::endl;
    double presamplings=60.0;
    double drifttick=(MCvertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198)+presamplings;
    
    std::cout<<"%%%%%%%%%%%%%%%%%%   drifttick= "<<std::setprecision(10)<<drifttick<<std::endl;
    ftime_vertex.push_back(drifttick);
    ftime_vertex.push_back(drifttick);
    
   

  }






}



  //////////////////////////////////////////////////////
  // here is how to get a collection of objects out of the file
  // and connect it to an art::Handle
  //////////////////////////////////////////////////////
  // Read in the clusterList object(s).
   
  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
  art::PtrVector<recob::EndPoint2D> endpointlist;
 
    for (unsigned int i = 0; i < endpointListHandle->size(); ++i){
      art::Ptr<recob::EndPoint2D> endpointHolder(endpointListHandle,i);
      endpointlist.push_back(endpointHolder);
    }
    
  std::cout<<"SHOULD BE GETTING RECO VERTEX, endpointlist.size()= "<<endpointlist.size()<<std::endl;
  
  if(endpointlist.size()==0){
  std::cout<<"ATTENTION: NO VERTEX FOUND, KINGACLUSTER WILL EXIT"<<std::endl;
  return;
  }
  for (unsigned int j = 0; j<endpointlist.size();j++){

     ftime_vertex_reco.push_back(endpointlist[j]->DriftTime());
     fwire_vertex_reco.push_back(endpointlist[j]->WireNum());
          
          std::cout<<"j="<<j<<" vtx2d_w="<<endpointlist[j]->WireNum()<<" vtx2d_t="<<endpointlist[j]->DriftTime()<<std::endl;
          
        }
      
     
 //********************************************************************     
       

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  //Point to a collection of clusters to output.
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);


  art::ServiceHandle<geo::Geometry> geom;
  art::PtrVector<recob::Cluster> clusIn;
 
  art::PtrVector<recob::Hit> hits;
  art::PtrVector<recob::Hit> clusterHits;
 
 
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
std::cout<<"No of DBSCAN clusters= "<<clusIn.size()<<std::endl;



 unsigned int p(0),w(0),t(0), channel(0);
 

 for(size_t tpc = 0; tpc < geom->NTPC(); ++tpc){
   std::cout<<"No of planes = "<<geom->Nplanes(tpc)<<std::endl;

   for(unsigned int plane = 0; plane < geom->Nplanes(tpc); plane++) {

    need_to_reassign_hitsIDs=0;
    go_ahead_at_reassign=0;
     for(unsigned int j=0; j<clusIn.size();++j) {
   
       hits=clusIn[j]->Hits();
       for(unsigned int i = 0; i< hits.size(); ++i){
	 channel=hits[i]->Wire()->RawDigit()->Channel();
	 geom->ChannelToWire(channel,t,p,w);

	 if(p == plane && t == tpc){
	   allhits.push_back(hits[i]);
	   //std::cout<<"plane= "<<plane<<" wire= "<<w<<" time= "<<hits[i]->PeakTime()<<std::endl; 
	 }
   
       }
       // std::cout<<"hits.size()= "<<hits.size()<<std::endl;
     }
   
     //std::cout<<"allhits.size()="<<allhits.size()<<std::endl;
   
   
     //Now we have hits for the plane that we are on right now, so let's do some work:
     //maxBin.clear();
     std::cout<<"ATTENTION, STARTING WORK ON PLANE# "<<plane<<std::endl;
     AngularDistribution(tpc,plane);
     FindMax(tpc,plane);
     if(fpeaks_found==0){
       std::cout<<"KingaClusters FAILED on this event because no peaks were found. Perhaps your threshold for peak's height is too big. Goodbye! "<<std::endl;
       allhits.clear();
       maxBin.clear();
       maxBinValues.clear();
       SortedMaxBin.clear();
       MaxStartPoint.clear();
       MaxEndPoint.clear();
       MaxStartPointTheta.clear();
       MaxEndPointTheta.clear();
       HitsWithClusterID.clear();
       FinalPeaks.clear();
       OriginalmaxBinValues.clear();
     for(int bin=0; bin< fh_theta_ind_2D->GetNbinsX(); bin++){
	 
	 fh_theta_ind_2D->SetBinContent(bin,0);
	 fh_theta_coll_2D->SetBinContent(bin,0);
	 fh_theta_ind->SetBinContent(bin,0);
	 fh_theta_coll->SetBinContent(bin,0);
	 fh_theta_coll_Area->SetBinContent(bin,0);
	 fh_theta_ind_Area->SetBinContent(bin,0);
       }
       
       return;
     }
     //FinalPeaks();
     FindClusters(tpc,plane);
     if(need_to_reassign_hitsIDs==1){
     std::cout<<"***************************************************************"<<std::endl;
     std::cout<<"***************  ATTENTION   ***********************"<<std::endl;
     std::cout<<" WILL NEED TO REASSIGN HIT IDs"<<std::endl;
     std::cout<<"***************************************************************"<<std::endl;
     FindClusters(tpc,plane);}
     std::cout<<"HitsWithClusterID.size()= "<<HitsWithClusterID.size()
	      << "compare with allhits.size()= "<<allhits.size()<<std::endl;
	      
	            
 //********************************************************************         
	      
	
	      
     for(unsigned int ClusterNo=0; ClusterNo<MaxStartPoint.size();ClusterNo++) {
    
       for(unsigned int j=0; j<HitsWithClusterID.size();j++){
     
       if(HitsWithClusterID[j]==(ClusterNo+1)){
       
       clusterHits.push_back(allhits[j]);
       } //if
    
    
    
    } //loop over HitsWithClusterID

// let's look at the clusters produced:
std::cout<<"For Cluster # "<<ClusterNo<<" we have "<<clusterHits.size()<<" hits :"<<std::endl;




//     //.................................
    if (clusterHits.size()>0)
	    {
	    
	    
	    
	    
	      /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	      unsigned int p = 0; 
	      unsigned int t = 0; 
	      unsigned int sw = 0;
	      unsigned int ew = 0;
	      geom->ChannelToWire(clusterHits[0]->Wire()->RawDigit()->Channel(), t, p, sw);
	      geom->ChannelToWire(clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel(), t, p, ew);

	      clusterHits.sort(cluster::SortByWire());
	      
	      recob::Cluster cluster(clusterHits, 
				     sw*1., 0.,
				     clusterHits[0]->PeakTime(), clusterHits[0]->SigmaPeakTime(),
				     ew*1., 0.,
				     clusterHits[clusterHits.size()-1]->PeakTime(), clusterHits[clusterHits.size()-1]->SigmaPeakTime(),
				     -999., 0., 
				     -999., 0.,
				    ClusterNo);
				    
				    
							    
std::cout<<"Produced Cluster #"<<ClusterNo<<std::endl;
	      ccol->push_back(cluster);
	      //std::cout<<"no of hits for this cluster is "<<clusterHits.size()<<std::endl;
	     // std::cout<<cluster.StartPos()[0]<<", "<<cluster.StartPos()[1]<<" --> "<<cluster.EndPos()[0]<<", "<<cluster.EndPos()[1]<<std::endl;
	      

	      
	      clusterHits.clear();
	      //////
	    }
      
   } //clusters
    
     allhits.clear();
     maxBin.clear();
     maxBinValues.clear();
     SortedMaxBin.clear();
     MaxStartPoint.clear();
     MaxEndPoint.clear();
     MaxStartPointTheta.clear();
     MaxEndPointTheta.clear();
     HitsWithClusterID.clear();
     FinalPeaks.clear();
     OriginalmaxBinValues.clear();
     std::cout<<"Should be starting to work on the other plane now"<<std::endl;
    }//Planes
    
 }// end loop over tpcs 
 
 
 evt.put(ccol);
 
 
   
 return;
}
    
//..............................................................  

 void cluster::KingaCluster::AngularDistribution(unsigned int tpc, unsigned int plane){   
 
 if(plane==0){
    for(int bin=0; bin< fh_theta_ind->GetNbinsX(); bin++){
   
    fh_theta_ind_2D->SetBinContent(bin,0);
    fh_theta_ind->SetBinContent(bin,0);
    fh_theta_ind_Area->SetBinContent(bin,0);
    }
   }
   
   if(plane==1){
    for(int bin=0; bin< fh_theta_ind_Area->GetNbinsX(); bin++){
   
    fh_theta_coll_2D->SetBinContent(bin,0);
    fh_theta_coll->SetBinContent(bin,0);
    fh_theta_coll_Area->SetBinContent(bin,0);
   
    }
   }
 
  
 art::ServiceHandle<geo::Geometry> geom;
 if(fMC==1){
   unsigned int channel2,plane2,wire2,tpc2;  
plane2=plane;

  if(plane==0){
	MCvertex[0]=.3;//force time coordinate to be closer to induction plane 
	}
      else{
	MCvertex[0]=-.3;//force time coordinate to be closer to collection plane
     }
      channel2 = geom->NearestChannel(MCvertex);
      geom->ChannelToWire(channel2,tpc2,plane2,wire2);   
   std::cout<<"%%%%%%%%%%%%%%%%%%   WIRE VERTEX IS: "<<wire2<<std::endl;
   fwire_vertex.push_back(wire2);
   
  
   
  } 
   
//std::vector<unsigned int> fwire_vertex,ftime_vertex;
double a_polar, b_polar,theta_polar;
 ftimetick      =  0.198; //get from parameterset
 fdriftvelocity =  0.157;  //get from paramtereset 9either k and V)
 fpi=3.141592653;


unsigned int channel=0, w=0;
 unsigned int p=plane, t=tpc;

//std::cout<<"for PLANE "<<plane<<" fwire_vertex= "<<fwire_vertex[plane]<<" ftime_vertex= "<<ftime_vertex[plane]<<std::endl;


std::cout<<"No of HITS for plane "<<plane<<" is: "<<allhits.size()<<std::endl;
 
  for(unsigned int i = 0; i< allhits.size(); ++i){
  
  // if(i==0){
//         fwire_vertex=allhits[i]->Wire()->RawDigit()->Channel();
//         ftime_vertex=allhits[i]->PeakTime();
//         std::cout<<"for PLANE "<<plane<<" fwire_vertex= "<<fwire_vertex<<" ftime_vertex= "<<ftime_vertex<<std::endl;
//           }   
 
 channel=allhits[i]->Wire()->RawDigit()->Channel();
 geom->ChannelToWire(channel,t,p,w);
 //std::cout<<"................................"<<std::endl;
 //std::cout<<" For hit# "<<i<<" w= "<<w<<" fwire_vertex[plane]= "<<fwire_vertex[plane]<<" ftime_vertex[plane]= "<<ftime_vertex[plane];
  
     int diff_w= w - fwire_vertex_reco[plane];
      b_polar = diff_w*0.4; /**in cm*/
      
      //std::cout<<" diff_w= "<<diff_w<<std::endl;
      //std::cout<<" b_polar= "<<b_polar<<std::endl;
      a_polar = (allhits[i]->PeakTime() - ftime_vertex_reco[plane])* ftimetick *fdriftvelocity; /** in cm*/
      
       
      theta_polar =fabs(asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2)))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/
      
     
     // We have 4 cases depending on which quater a hit is (origin being defined on a vertex):
     
     if(b_polar==0 && a_polar==0){
     theta_polar = 90;/** in deg*/
      }
     else if(b_polar==0 && a_polar<0){
     theta_polar = 180;/** in deg*/
      }
     else if(b_polar==0 && a_polar>0){
     theta_polar = 0;/** in deg*/
      }
     else if(b_polar>0 && a_polar==0){
     theta_polar = 90;/** in deg*/
      }
      else if(b_polar<0 && a_polar==0){
     theta_polar = -90;/** in deg*/
      }
      
      
     
     else if(b_polar>0 && a_polar>0){
     theta_polar = 90-theta_polar;/** in deg*/
       /** in deg*/
      }
      else if(b_polar>0 && a_polar<0){
       theta_polar = 90+theta_polar;/** in deg*/
      }
      else if(b_polar<0 && a_polar>0){
       theta_polar = -(90-theta_polar);/** in deg*/
      }
      else if(b_polar<0 && a_polar<0){
       theta_polar = -(90+theta_polar);/** in deg*/
      }
      
     //fh_theta[plane]->Fill(theta_polar,allhits[i]->Charge());
     //std::cout<<"**** theta_polar= "<<theta_polar<<std::endl;
    // std::cout<<" a_polar= "<<a_polar<<" b_polar= "<<b_polar<<std::endl;
     
     //std::cout<<"w= "<<w<<" t= "<<std::setprecision(10)<<allhits[i]->PeakTime();
     //std::cout<<" theta_polar= "<<theta_polar<<" hit area= "<<(allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4<<std::endl;
if (plane==0 ) {

//std::cout<<"plane= "<<plane<<" theta= "<<theta_polar<<"  channel= "<<allhits[i]->Wire()->RawDigit()->Channel()<<" time= "<<allhits[i]->PeakTime()<<" fwire_vertex[plane]= "<<fwire_vertex[plane]<<" ftime_vertex[plane]= "<<ftime_vertex[plane]<<std::endl;


fh_theta_ind->Fill(theta_polar);
fh_theta_ind_2D->Fill(theta_polar,allhits[i]->Charge());
 fh_theta_ind_Area->Fill(theta_polar,(allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);
 
 Hit_Area_Ind->Fill((allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);
}
if (plane==1 ) {

//std::cout<<"plane= "<<plane<<"  w= "<<w<<" time= "<<allhits[i]->PeakTime()<<" theta= "<<theta_polar<<" fwire_vertex[plane]= "<<fwire_vertex[plane]<<" ftime_vertex[plane]= "<<ftime_vertex[plane];


// if(w>95 && w<118 && allhits[i]->PeakTime()>500 && allhits[i]->PeakTime()<800){
// std::cout<<"  ***********"<<std::endl;
// }
// else{ std::cout<<std::endl;}




fh_theta_coll->Fill(theta_polar);
fh_theta_coll_2D->Fill(theta_polar,allhits[i]->Charge());
fh_theta_coll_Area->Fill(theta_polar,(allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);

Hit_Area_Coll->Fill((allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);

}
  }
 
  // for(int i=0;i<startTimes.size();i++){
//   std::cout<<"startTime is at bin = "<<startTimes[i]<<" and its value is "<<fh_theta_ind->GetBinContent(startTimes[i])<<std::endl;
//   }
//   for(int i=0;i<endTimes.size();i++){
//   std::cout<<"endTime is at bin = "<<endTimes[i]<<" and its value is "<<fh_theta_ind->GetBinContent(endTimes[i])<<std::endl;
//   }
 
 
 
    }
    
    
    
//..............................................................   
 void cluster::KingaCluster::FindMax(unsigned int tpc, unsigned int plane){  

 // std::cout<<"No of bins= "<<fh_theta_ind->GetNbinsX()<<std::endl;
//   std::cout<<" Bincontent= "<<fh_theta_ind->GetBinContent(48)<<std::endl;
//   std::cout<<" Bincontent= "<<fh_theta_ind->GetBinContent(49)<<std::endl;
//   std::cout<<" Bincontent= "<<fh_theta_ind->GetBinContent(50)<<std::endl;
// std::vector<int> PossibleFinalMax, FinalMax;
   std::vector<int> startTimes;  //stores time of 1st local minimum
    // std::vector<int> maxBin;    //stores time of local maximum
    std::vector<int> endTimes;    //stores time of 2nd local minimum
   bool maxFound; //Flag for whether a value>threshold has been found
  
   int minTimeHolder;
 
  startTimes.clear();
  maxBin.clear();
  endTimes.clear();
 //std::cout<<"We have "<<allhits.size()<<" hits for plane "<<plane<<std::endl;

 //  double threshold=76*allhits.size()/2;
//   double MinThreshold=1000;
//  double threshold=6;
//   double MinThreshold=4;

// worked for pizero and nu_e :
 // double threshold=80000;
//   double MinThreshold=60000;
  //............................
  
  
  
  // double threshold=30000;
//   double MinThreshold=20000;

  //for lines:
  // double threshold=20000;
//   double MinThreshold=10000;

  // double threshold=600;
//   double MinThreshold=100;

  double threshold=0.2;
  //double MinHitThreshold=1; //Making sure that the peak found in fh_theta_coll_Area and _ind corresponds to more than 1 hit. Sometimes you can have a very large hit area that will produce a peak in that distribution but it only corresponds to 1 hit.
  double MinThreshold=0.1;
  
  
  int time=1;
  int ValidPeak=0;
  std::cout<<"Threshold that a peak must have in order to be considered a peak = "<<threshold<<". For inflection points we must have the value to drop to "<<MinThreshold<<std::endl;
  
 //......................................................... 
 //collection plane:
 if (plane==1){
 

  for(int bin=1; bin<fh_theta_coll_Area->GetNbinsX()+1;bin++){
 
  if(fh_theta_coll_Area->GetBinContent(bin)>fh_theta_coll_Area->GetBinContent(bin+1) && fh_theta_coll_Area->GetBinContent(bin+1)<fh_theta_coll_Area->GetBinContent(bin+2)) {
      //only add points if we've already found a local max above threshold.
      if(maxFound) {
        endTimes.push_back(time+1);
        maxFound = false;
        //keep these in case new hit starts right away
        minTimeHolder = time+2;
       }
      else {minTimeHolder = time+1; }
    }
  //if not a minimum,-> test if we are at a local maximum
    //if so and the max value is above threshold add it and proceed.
    else if(fh_theta_coll_Area->GetBinContent(bin)<fh_theta_coll_Area->GetBinContent(bin+1) &&
        fh_theta_coll_Area->GetBinContent(bin+1)>fh_theta_coll_Area->GetBinContent(bin+2) &&
        fh_theta_coll_Area->GetBinContent(bin+1) > threshold ) {
      maxFound = true;
      maxBin.push_back(time+1);
      startTimes.push_back(minTimeHolder);         
    }

    time++; 
     
  }
  
  if(maxBin.size()==0){
  
  std::cout<<"COLLECTION PLANE: "<<std::endl;
  std::cout<<" COULDN'T FIND ANY MAXIMA IN YOUR THETA DISTRIBUTION!!!! PROGRAM WILL NOT PROCEED!!!"<<std::endl;
    fpeaks_found=0;}
 if(maxBin.size()>0){     
     
 // for(unsigned int i=0;i<maxBin.size();i++){
//   std::cout<<"maxTime is at bin = "<<maxBin[i]<<" and its value is "<<fh_theta_coll_Area->GetBinContent(maxBin[i])<<std::endl;
//   }
  
  
  // Lets make sure that the first bin in the maxBin corresponds to the highest peak:

//std::vector<double> maxBinValues;
  for(unsigned int i=0;i<maxBin.size();i++){
  maxBinValues.push_back(fh_theta_coll_Area->GetBinContent(maxBin[i]));
  OriginalmaxBinValues.push_back(fh_theta_coll_Area->GetBinContent(maxBin[i]));
  }
  //std::cout<<"The largest is at position:  "<<std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))<<" which corresponds to bin #   "<<maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]<<" and its value= "<<*std::max_element(maxBinValues.begin(),maxBinValues.end())<<std::endl;
  
  //sort values from the largest to the smallest in maxBinValues, then find the corresponding bin numbers to create SortedMaxBin:
  
  sort(maxBinValues.begin(),maxBinValues.end());
//std::cout<<"maxBinValues after sort:"<<std::endl;
 // for(unsigned int i=0;i<maxBinValues.size();i++){
// 
//    std::cout<<maxBinValues[i]<<std::endl;
//  }

reverse (maxBinValues.begin(),maxBinValues.end());
 //std::cout<<"maxBinValues in the correct order are now:"<<std::endl;
 // for(unsigned int i=0;i<maxBinValues.size();i++){
// 
//    std::cout<<maxBinValues[i]<<std::endl;
//  }

 for(unsigned int i=0; i<maxBinValues.size();i++){

   std::vector<double>::iterator pos=std::find( OriginalmaxBinValues.begin(), OriginalmaxBinValues.end(),maxBinValues[i]);
   SortedMaxBin.push_back(maxBin[pos-OriginalmaxBinValues.begin()]);

 }


  //create SortedMaxBin vector whose first element is the global max:
  
 //  SortedMaxBin.push_back(maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]);
//   for(int i=0; i<maxBin.size();i++){
//    if(i!=std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end())))
//   SortedMaxBin.push_back(maxBin[i]);
  
//   }
  
 //  std::cout<<"SortexMaxBin elements are: "<<std::endl;
// for(unsigned int i=0; i<SortedMaxBin.size(); i++)
// {
// 
// std::cout<<SortedMaxBin[i]<<std::endl;
// 
// }
// int ValidPeak=0;
  // loop over maxima and find where they start on the left and right side, form cluster for each:
  
  for(unsigned int maxNo=0; maxNo<SortedMaxBin.size();maxNo++){
  
  //loop over the ranges and make sure that your peaks don't fall into already formed clusters
  //std::cout<<"Right now we have "<<MaxStartPoint.size()<<" ranges"<<std::endl;
  if(MaxStartPoint.size()==0){ValidPeak=1;}
  for(unsigned int NoRange=0; NoRange<MaxStartPoint.size();NoRange++){
  //std::cout<<"Checking peak "<<SortedMaxBin[maxNo]<<std::endl;
  if(SortedMaxBin[maxNo]>MaxStartPoint[NoRange] && SortedMaxBin[maxNo]<MaxEndPoint[NoRange]){
  //maxNo++;
  //std::cout<<"this peak is out of the picture! --> "<<SortedMaxBin[maxNo]<<std::endl;
  ValidPeak=0;
  break;}
  else {ValidPeak=1;
  //std::cout<<"passed"<<std::endl;
  }
  }
   if(ValidPeak==1){
   
  // std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
  FinalPeaks.push_back(SortedMaxBin[maxNo]);
   std::cout<<"We are working on peak at bin #"<<SortedMaxBin[maxNo]<<std::endl;
  //start at the peak and go left
      for(int LeftBin=SortedMaxBin[maxNo]-1;LeftBin>SortedMaxBin[maxNo]-30; LeftBin--)
      {
        if(fh_theta_coll_Area->GetBinContent(LeftBin)<MinThreshold && fh_theta_coll_Area->GetBinContent(LeftBin-1)>fh_theta_coll_Area->GetBinContent(LeftBin)) {
           MaxStartPoint.push_back(LeftBin);
          // std::cout<<"picked option 1, startin point @bin "<<LeftBin<<std::endl;
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin<<"("<<-180+2*LeftBin<<" degrees)"<<" RightBin= ";
           break;
           }
           else if(fh_theta_coll_Area->GetBinContent(LeftBin)<MinThreshold && fh_theta_coll_Area->GetBinContent(LeftBin-1)==0){
           MaxStartPoint.push_back(LeftBin-1);
           //std::cout<<"picked option 2, startin point @bin "<<LeftBin-1<<std::endl;
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin-1<<"("<<-180+2*(LeftBin-1)<<" degrees)"<<" RightBin= ";
           break;
           
           }
           else if (LeftBin==SortedMaxBin[maxNo]-29){std::cout<<" cannot find starting point of the peak!!!!!"<<std::endl;}
         
      }
     
   
   
   
  
  for(int RightBin=SortedMaxBin[maxNo]+1;RightBin<SortedMaxBin[maxNo]+30; RightBin++)
      {
        if(fh_theta_coll_Area->GetBinContent(RightBin)<MinThreshold && fh_theta_coll_Area->GetBinContent(RightBin+1)>fh_theta_coll_Area->GetBinContent(RightBin)) {
           MaxEndPoint.push_back(RightBin);
           std::cout<<RightBin<<"("<<-180+2*RightBin<<" degrees)"<<std::endl;
           break;
           }
         else if(fh_theta_coll_Area->GetBinContent(RightBin)<MinThreshold && fh_theta_coll_Area->GetBinContent(RightBin+1)==0){
           MaxEndPoint.push_back(RightBin+1);
           std::cout<<RightBin+1<<"("<<-180+2*(RightBin+1)<<" degrees)"<<std::endl;
           break;
           }
         else if(RightBin==SortedMaxBin[maxNo]+29){std::cout<<" cannot find end point of the peak!!!!!"<<std::endl;}
         
      }
      
      
      
      
  } //valid peak
  
  ValidPeak=0;
  
}//peaks
  
  

  } //if maxBin.size()>0
  
  startTimes.clear();
  maxBin.clear();
  endTimes.clear();
  
 } //plane 1
 
  time=1;
  ValidPeak=0;
  
 //......................................................... 
 
 
 //induction plane
 if(plane==0){
 
 
 //std::cout<<"No of bins= "<<fh_theta_ind_Area->GetNbinsX()<<std::endl;
for(int bin=1; bin<fh_theta_ind_Area->GetNbinsX()+1;bin++){
 
  if(fh_theta_ind_Area->GetBinContent(bin)>fh_theta_ind_Area->GetBinContent(bin+1) && fh_theta_ind_Area->GetBinContent(bin+1)<fh_theta_ind_Area->GetBinContent(bin+2)) {
      //only add points if we've already found a local max above threshold.
      if(maxFound) {
        endTimes.push_back(time+1);
        maxFound = false;
        //keep these in case new hit starts right away
        minTimeHolder = time+2;
       }
      else {minTimeHolder = time+1; }
    }
  //if not a minimum test if we are at a local maximum
    //if so and the max value is above threshold add it and proceed.
    else if(fh_theta_ind_Area->GetBinContent(bin)<fh_theta_ind_Area->GetBinContent(bin+1) &&
        fh_theta_ind_Area->GetBinContent(bin+1)>fh_theta_ind_Area->GetBinContent(bin+2) &&
        fh_theta_ind_Area->GetBinContent(bin+1) > threshold) {
      maxFound = true;
      maxBin.push_back(time+1);
      startTimes.push_back(minTimeHolder); 
      
      // && fh_theta_ind->GetBinContent(bin+1) > MinHitThreshold
    }
    time++;
 
  }
 std::cout<<"INDUCTION PLANE: "<<std::endl;

  if(maxBin.size()==0){
  std::cout<<" COULDN'T FIND ANY MAXIMA IN YOUR THETA DISTRIBUTION!!!! PROGRAM WILL NOT PROCEED!!!"<<std::endl;
    fpeaks_found=0;}
    
 if(maxBin.size()>0){   
 std::cout<<"maxBin.size()= "<<maxBin.size()<<std::endl;
  for(unsigned int i=0;i<maxBin.size();i++){
  std::cout<<"maxTime is at bin = "<<maxBin[i]<<" ("<<-180+2*maxBin[i]<<" degrees)"<<" and its value is "<<fh_theta_ind_Area->GetBinContent(maxBin[i])<<std::endl;
std::cout<<"...................................."<<std::endl;
 }//maxBin
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////





// Lets make sure that the first bin in the maxBin corresponds to the highest peak:

//std::vector<double> maxBinValues;
  for(unsigned int i=0;i<maxBin.size();i++){
  maxBinValues.push_back(fh_theta_ind_Area->GetBinContent(maxBin[i]));
  OriginalmaxBinValues.push_back(fh_theta_ind_Area->GetBinContent(maxBin[i]));

  }
  //std::cout<<"The largest is at position:  "<<std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))<<" which corresponds to bin #   "<<maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]<<" and its value= "<<*std::max_element(maxBinValues.begin(),maxBinValues.end())<<std::endl;
  
  //sort values from the largest to the smallest in maxBinValues, then find the corresponding bin numbers to create SortedMaxBin:
  
  sort(maxBinValues.begin(),maxBinValues.end());
 // std::cout<<"maxBinValues after sort:"<<std::endl;
//  for(unsigned int i=0;i<maxBinValues.size();i++){
// 
//    std::cout<<maxBinValues[i]<<std::endl;
//  }

reverse (maxBinValues.begin(),maxBinValues.end());
 // std::cout<<"maxBinValues in the correct order are now:"<<std::endl;
//  for(unsigned int i=0;i<maxBinValues.size();i++){
// 
//    std::cout<<maxBinValues[i]<<std::endl;
//  }

 for(unsigned int i=0; i<maxBinValues.size();i++){

   std::vector<double>::iterator pos=std::find( OriginalmaxBinValues.begin(), OriginalmaxBinValues.end(),maxBinValues[i]);
   SortedMaxBin.push_back(maxBin[pos-OriginalmaxBinValues.begin()]);

 }
  
  
  
  
  
  
  
  
  //create SortedMaxBin vector whose first element is the global max:
  
  
  // std::vector<int> SortedMaxBin;
  // SortedMaxBin.push_back(maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]);
//   for(int i=0; i<maxBin.size();i++){
//    if(i!=std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end())))
//   SortedMaxBin.push_back(maxBin[i]);
//   
//   }
  
 //  std::cout<<"SortexMaxBin elements are: "<<std::endl;
// for(unsigned int i=0; i<SortedMaxBin.size(); i++)
// {
// 
// std::cout<<SortedMaxBin[i]<<std::endl;
// 
// }

  // loop over maxima and find where they start on the left and right side, form cluster for each:
  
  for(unsigned int maxNo=0; maxNo<SortedMaxBin.size();maxNo++){
  
  //loop over the ranges and make sure that your peaks don't fall into already formed clusters
  if(MaxStartPoint.size()==0){ValidPeak=1;}
  for(unsigned int NoRange=0; NoRange<MaxStartPoint.size();NoRange++){
  if(SortedMaxBin[maxNo]>MaxStartPoint[NoRange] && SortedMaxBin[maxNo]<MaxEndPoint[NoRange]){
  ValidPeak=0;
  break;}
  else {ValidPeak=1;}
  }
  
   if(ValidPeak==1){
   
  //std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
  FinalPeaks.push_back(SortedMaxBin[maxNo]);
  //start at the peak and go left
      for(int LeftBin=SortedMaxBin[maxNo]-1;LeftBin>SortedMaxBin[maxNo]-30; LeftBin--)
      {
        if(fh_theta_ind_Area->GetBinContent(LeftBin)<MinThreshold && fh_theta_ind_Area->GetBinContent(LeftBin-1)>fh_theta_ind_Area->GetBinContent(LeftBin)) {
           MaxStartPoint.push_back(LeftBin);
            
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin<<"("<<-180+2*LeftBin<<" degrees)"<<" RightBin= ";
        
           break;
           }
           else if(fh_theta_ind_Area->GetBinContent(LeftBin)<MinThreshold && fh_theta_ind_Area->GetBinContent(LeftBin-1)==0){
           MaxStartPoint.push_back(LeftBin-1);
           //std::cout<<"picked option 2, startin point @bin "<<LeftBin-1<<std::endl;
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin-1<<"("<<-180+2*(LeftBin-1)<<" degrees)"<<" RightBin= ";
           break;
           
           }
           else if (LeftBin==SortedMaxBin[maxNo]-30){std::cout<<"cannot find starting point of the peak!!!!!"<<std::endl;}
         
      }
  
  
  for(int RightBin=SortedMaxBin[maxNo]+1;RightBin<SortedMaxBin[maxNo]+30; RightBin++)
      {
        if(fh_theta_ind_Area->GetBinContent(RightBin)<MinThreshold && fh_theta_ind_Area->GetBinContent(RightBin+1)>fh_theta_ind_Area->GetBinContent(RightBin)) {
           MaxEndPoint.push_back(RightBin);
           std::cout<<RightBin<<"("<<-180+2*RightBin<<" degrees)"<<std::endl;
           break;
           }
           else if(fh_theta_ind_Area->GetBinContent(RightBin)<MinThreshold && fh_theta_ind_Area->GetBinContent(RightBin+1)==0){
           MaxEndPoint.push_back(RightBin+1);
           std::cout<<RightBin+1<<"("<<-180+2*(RightBin+1)<<" degrees)"<<std::endl;
           break;
           }
           else if(RightBin==SortedMaxBin[maxNo]+20){std::cout<<"cannot find end point of the peak!!!!!"<<std::endl;}
         
      }
  
  }
  ValidPeak=0;
  
}//peaks
  
  
}// if maxBin.size()>0 this means that we can find peaks in the theta distribution 
  
  
 } //plane 0
 //......................................................... 
 //Now clear the histograms:
 
 // for(int bin=0; bin< fh_theta_ind_2D->GetNbinsX(); bin++){
//    
//    fh_theta_ind_2D->SetBinContent(bin,0);
//    fh_theta_coll_2D->SetBinContent(bin,0);
//    fh_theta_ind->SetBinContent(bin,0);
//    fh_theta_coll->SetBinContent(bin,0);
//    fh_theta_coll_Area->SetBinContent(bin,0);
//    fh_theta_ind_Area->SetBinContent(bin,0);
//  }
 
}   

//..............................................................   
// void cluster::KingaCluster::FinalPeaks(){  
// std::cout<<"In FinalPeaks()"<<std::endl;
// // for(int i=0;i<maxBin.size();i++){
// //   std::cout<<"maxTime is at bin = "<<maxBin[i]<<" and its value is "<<fh_theta_ind_Area->GetBinContent(maxBin[i])<<std::endl;
// //   }
// 
// 
// 
// }

//..............................................................   


void cluster::KingaCluster::FindClusters(unsigned int tpc, unsigned int plane){ 

//First order check: make sure ranges of the peaks don't overlap, if they do you need to correct for that. The highest peaks should be left alone and we should work with the ones of the smallest peak value first. This means start from the end of the sorted peaks.

std::vector<int> peak_with_wrong_range, peak_with_which_it_ovelaps;
peak_with_wrong_range.clear();
peak_with_which_it_ovelaps.clear();



for( int pk=FinalPeaks.size()-1; pk>=0;pk--){
   std::cout<<"pk= "<<pk<<std::endl;
   for(unsigned int pk2=0; pk2<FinalPeaks.size();pk2++){
    if(pk!= (int) pk2 && ((MaxStartPoint[pk]<MaxEndPoint[pk2] && MaxStartPoint[pk]>MaxStartPoint[pk2])||( MaxEndPoint[pk]>MaxStartPoint[pk2] && MaxEndPoint[pk]<MaxEndPoint[pk2] ))){
    std::cout<<"WRONG RANGE, NEED TO FIX IT FOR PEAK AT BIN #"<<FinalPeaks[pk]<<std::endl;
    
    peak_with_wrong_range.push_back(pk); //this gives peak#, NOT a bin#
    peak_with_which_it_ovelaps.push_back(pk2); //this gives peak#, NOT a bin#
    }
  
   }


}

int diff=0;


for(unsigned int i=0; i<peak_with_wrong_range.size(); i++){


//if front of a range overlaps with the back of the range already in place
if(MaxStartPoint[peak_with_wrong_range[i]]<MaxEndPoint[peak_with_which_it_ovelaps[i]] && MaxStartPoint[peak_with_wrong_range[i]]>MaxStartPoint[peak_with_which_it_ovelaps[i]]){
 
 diff=MaxEndPoint[peak_with_which_it_ovelaps[i]]-MaxStartPoint[peak_with_wrong_range[i]];
  MaxStartPoint[peak_with_wrong_range[i]]=MaxStartPoint[peak_with_wrong_range[i]]+ diff;
  std::cout<<"changing startpoint to bin #"<<MaxStartPoint[peak_with_wrong_range[i]]<<std::endl;

}

//if back of a range overlaps with the front of the range already in place
if(MaxEndPoint[peak_with_wrong_range[i]]>MaxStartPoint[peak_with_which_it_ovelaps[i]] && MaxEndPoint[peak_with_wrong_range[i]]<MaxEndPoint[peak_with_which_it_ovelaps[i]]){
 
 diff=MaxEndPoint[peak_with_wrong_range[i]]-MaxStartPoint[peak_with_which_it_ovelaps[i]];
  MaxEndPoint[peak_with_wrong_range[i]]=MaxEndPoint[peak_with_wrong_range[i]]- diff;
  
  std::cout<<"changing endpoint to bin #"<<MaxEndPoint[peak_with_wrong_range[i]]<<std::endl;


}



diff=0;


}






// First let's make sure that each range of peaks contains more than some specified number of hits. You can do it by knowing the ranges, going into histograms and counting the entries from start of the range to the end.If some range contains less than the specified number you need to remove the peak.

std::cout<<" NO OF FINALPEAKS BEFORE EVALUATION IS: "<<FinalPeaks.size()<<std::endl;





int MinHitsInRange=2; //later make it a parameter. DEFINED BELOW, LOOK->>>
double no_hits_in_range=0;
std::vector<int> TempFinalPeaks,TempMaxStartPoint,TempMaxEndPoint;
TempFinalPeaks.clear();
TempMaxStartPoint.clear();
TempMaxEndPoint.clear();

std::vector<int> positive_diff_end_minus_start;
std::vector<int> positive_diff_start_minus_end;
std::vector<int> diff_end_minus_start;
std::vector<int> diff_start_minus_end;

 positive_diff_end_minus_start.clear();
 positive_diff_start_minus_end.clear();
 diff_end_minus_start.clear();
 diff_start_minus_end.clear();
double closest_range_right_side=0;
double closest_range_left_side=0;
int this_is_the_first_range=0;
int this_is_the_last_range=0;
int well_separated=0;

//---------------START BASIC EVALUATION FIRST---------------




//Now I will eliminate all the peaks that have only 1 hit in it's RANGE(not in actual peak):
//However, if there is a group of peaks that have just one hit in their range and are all close to each other then we need to form a cluster out of them, call it "hand_made_peak". But this peak has to be well seprated from the rest!!!! so need to determine it's range also in order to check the separation

std::vector<int> one_hit_peaks,one_hit_peaks_start_point,one_hit_peaks_end_point,double_hit_peaks,double_hit_peaks_start_point,double_hit_peaks_end_point;
one_hit_peaks.clear();
one_hit_peaks_start_point.clear();
one_hit_peaks_end_point.clear();
double_hit_peaks.clear();
double_hit_peaks_start_point.clear();
double_hit_peaks_end_point.clear();

for(unsigned int peak=0; peak<MaxStartPoint.size(); peak++){
  for(int bin=MaxStartPoint[peak]; bin<MaxEndPoint[peak];bin++){
 
 
 //std::cout<<"bin= "<<bin<<std::endl;
   if(plane==0){no_hits_in_range+=fh_theta_ind->GetBinContent(bin);
 //std::cout<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range<<std::endl;
   }
   if(plane==1){no_hits_in_range+=fh_theta_coll->GetBinContent(bin);
 //std::cout<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range<<std::endl;
   }
 
  }//loop thru bins for each peak
  
  if(no_hits_in_range==1){
  one_hit_peaks.push_back(FinalPeaks[peak]);
  one_hit_peaks_start_point.push_back(MaxStartPoint[peak]);
  one_hit_peaks_end_point.push_back(MaxEndPoint[peak]);
  }
  
  if(no_hits_in_range==2){
  double_hit_peaks.push_back(FinalPeaks[peak]);
  double_hit_peaks_start_point.push_back(MaxStartPoint[peak]);
  double_hit_peaks_end_point.push_back(MaxEndPoint[peak]);
  }
  
  
  std::cout<<"no_hits_in_range= "<<no_hits_in_range<<" for peak at bin# "<<FinalPeaks[peak]<<std::endl;
  if(no_hits_in_range>2){
  TempFinalPeaks.push_back(FinalPeaks[peak]);
  TempMaxStartPoint.push_back(MaxStartPoint[peak]);
  TempMaxEndPoint.push_back(MaxEndPoint[peak]);
 }
  no_hits_in_range=0;
} //loop thru all peaks

FinalPeaks=TempFinalPeaks;
MaxStartPoint=TempMaxStartPoint;
MaxEndPoint=TempMaxEndPoint;

TempFinalPeaks.clear();
TempMaxStartPoint.clear();
TempMaxEndPoint.clear();
no_hits_in_range=0;

std::cout<<" So now the size of FinalPeaks.size()="<<FinalPeaks.size()<<std::endl;

//First, let's look at 2-hit peaks and a case in which there are TWO 2-hit peaks right next to each other. If this is the case make one peak out of them by increasing the range of the peak, make the peak be in the middle of the range.

//.......2-hit peaks START................................

// std::vector<int> used_double_hit_peaks;
// used_double_hit_peaks.clear();
// 
//   for(unsigned int i=0;i<double_hit_peaks.size(); i++){
//     for(unsigned int j=0;j<double_hit_peaks.size(); j++){
// 
//   if(i!=j && double_hit_peaks_end_point[i]==double_hit_peaks_start_point[j]){
//   
//   //create one peak with range of both peaks:
//   //Also make sure you do the combination only once (you can have i=1 and j=2 && //i=2 and j=1...you should just do one since they are the same)
//   
//   if(std::find(FinalPeaks.begin(),FinalPeaks.end(),(double_hit_peaks[i]+double_hit_peaks[j])*0.5)==FinalPeaks.end()){
//   FinalPeaks.push_back((double_hit_peaks[i]+double_hit_peaks[j])*0.5);
//   MaxStartPoint.push_back(double_hit_peaks_start_point[i]);
//   MaxEndPoint.push_back(double_hit_peaks_end_point[j]);
//   
//   used_double_hit_peaks.push_back(i);
//   used_double_hit_peaks.push_back(j);
//   std::cout<<"Creating one peak from TWO 2-hit peaks at peak at bin #"<<(double_hit_peaks[i]+double_hit_peaks[j])*0.5<<" and range in bins: ["<<double_hit_peaks_start_point[i]<<", "<<double_hit_peaks_end_point[j]<<"]"<<std::endl;
//   
//   }//if this peak doesn't exist yet, make it
//   
//   }
// 
// 
// 
// 
//     }
//    }

//Now fill FinalPeaks with the rest of double_hit_peaks, unless you can form a bigger cluster out of very close-by 2-hit-clusters:


//............. GROUP 2-HIT-CLUSTERS HERE, THE SAME WAY AS 1-HIT-CLUSTERS( see work on 1-hit clusters below)

// std::vector<int> temp_double_hit_peaks;
// temp_double_hit_peaks.clear();
// std::vector<int> temp_double_hit_peaks_start_point;
// temp_double_hit_peaks_start_point.clear();
// std::vector<int> temp_double_hit_peaks_end_point;
// temp_double_hit_peaks_end_point.clear();


// for(unsigned int k=0;k<double_hit_peaks.size();k++){
// 
// if(std::find(used_double_hit_peaks.begin(),used_double_hit_peaks.end(),double_hit_peaks[k])==used_double_hit_peaks.end()){
// 
// temp_double_hit_peaks.push_back(double_hit_peaks[k]);
// temp_double_hit_peaks_start_point.push_back(double_hit_peaks_start_point[k]);
// temp_double_hit_peaks_end_point.push_back(double_hit_peaks_end_point[k]);
// }
// }
std::vector<int> grouped_double_hit_peaks;
grouped_double_hit_peaks.clear();
// double_hit_peaks.clear();
// double_hit_peaks_start_point.clear();
// double_hit_peaks_end_point.clear();
// double_hit_peaks=temp_double_hit_peaks;
// double_hit_peaks_start_point=temp_double_hit_peaks_start_point;
// double_hit_peaks_end_point=temp_double_hit_peaks_end_point;
int diff_between_double_hit_peaks=7;
std::vector<int> used_peak;
used_peak.clear();
if(double_hit_peaks.size()>=2){

  for(unsigned int i=0;i<double_hit_peaks.size(); i++){
    for(unsigned int j=0;j<double_hit_peaks.size(); j++){
    
      if(grouped_double_hit_peaks.size()==0){
      
      // loop thru existing ranges and make sure that your 2-hit peaks are close to each other but also well separated from the existing ranges
       for(unsigned int r=0; r<FinalPeaks.size(); r++){ 
       std::cout<<"### FinalPeaks.size()= "<<FinalPeaks.size()<<std::endl;
       std::cout<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]<<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r]<<std::endl;
       if(i!=j && abs(double_hit_peaks[i]-double_hit_peaks[j])<=diff_between_double_hit_peaks ){ 
       
       //std::cout<<"abs(one_hit_peaks[i]-MaxStartPoint[r])= "<<abs(one_hit_peaks[i]-MaxStartPoint[r])<<" :: "<<one_hit_peaks[i]<<" - "<<MaxStartPoint[r]<<std::endl;
       
       
       //std::cout<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]<<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r]<<std::endl;
       
       //std::cout<<"one_hit_peaks[i]= "<<one_hit_peaks[i]<<" one_hit_peaks[j]= "<<one_hit_peaks[j]<<std::endl;
       
       
      //if you can't find this peak in the group, then add it:
         if(std::find(grouped_double_hit_peaks.begin(),grouped_double_hit_peaks.end(),double_hit_peaks[i])==grouped_double_hit_peaks.end()){
      grouped_double_hit_peaks.push_back(double_hit_peaks[i]);
      std::cout<<"ADDING 2-HIT PEAK #"<<double_hit_peaks[i]<<" TO THE GROUP"<<std::endl;
      
      used_peak.push_back(double_hit_peaks[i]);
         }
        if(std::find(grouped_double_hit_peaks.begin(),grouped_double_hit_peaks.end(),double_hit_peaks[j])==grouped_double_hit_peaks.end()){
     grouped_double_hit_peaks.push_back(double_hit_peaks[j]);
     std::cout<<"ADDING 2-HIT PEAK #"<<double_hit_peaks[j]<<" TO THE GROUP"<<std::endl;
     
     used_peak.push_back(double_hit_peaks[j]);
         }
   } //if difference between peak bin#s is small 
   }// loop thru already existing peaks 
} //if size=0 
    
     
    

     } //first loop thru each 1-hit peak
     
     
     
     
  } //second loop thru each 1-hit peak
  
  int failed=0;
   //Now go thru all the peaks again to find the third matching peak:
     if(grouped_double_hit_peaks.size()>0){
      for(unsigned int y=0;y<double_hit_peaks.size();y++){
       //now every added peak must be very close to EACH already added peak in the group
       for(unsigned int q=0; q<grouped_double_hit_peaks.size();q++){
       
        if(abs(double_hit_peaks[y]-grouped_double_hit_peaks[q])>5){
        failed=1;
        break;
        }
       
       } //for
       
       
       if(failed==0 && std::find(grouped_double_hit_peaks.begin(),grouped_double_hit_peaks.end(),double_hit_peaks[y])==grouped_double_hit_peaks.end() ){
       
       grouped_double_hit_peaks.push_back(double_hit_peaks[y]);
       
       used_peak.push_back(double_hit_peaks[y]);
       
       }
       failed=0;
       } //loop thru all 2-hit peaks
     }// if size >0
  
  
}


//so now we should have a group of 2-hit peaks that are close to each other. 

if(grouped_double_hit_peaks.size()>=2){


 //form a peak by taking the average of them and picking the one that's closest to that value:
 int sum=0;
 int made_peak=200;
 
 
   for(unsigned int k=0; k<grouped_double_hit_peaks.size(); k++){
    
    sum+=grouped_double_hit_peaks[k];
  // std::cout<<"grouped_double_hit_peaks.size()= "<<grouped_double_hit_peaks.size()<<" sum= "<<sum<<std::endl;
   }
   
    made_peak=sum/grouped_double_hit_peaks.size();
    std::cout<<"---------MADE A HOME_MADE_PEAK FROM GROUPS OF ***2***-HIT CLUSTERS--------THIS PEAK IS AT BIN #"<<made_peak;
    FinalPeaks.push_back(made_peak);
    //Now work on the range. For now just start at the smallest peak bin-1 and end at the largest peak bin+1:
    
    std::sort(grouped_double_hit_peaks.begin(),grouped_double_hit_peaks.end());
    
    MaxStartPoint.push_back(grouped_double_hit_peaks[0]-1);
    MaxEndPoint.push_back(grouped_double_hit_peaks[grouped_double_hit_peaks.size()-1]+1);
    
    std::cout<<" its range is ["<<grouped_double_hit_peaks[0]-1<<", "<<grouped_double_hit_peaks[grouped_double_hit_peaks.size()-1]+1<<"]"<<std::endl;
    
    
    
}


//If we weren't able to group 2-hit-clusters together then just add them to FinalPeaks the way they are:

for(unsigned int f=0; f<double_hit_peaks.size();f++){

 if(std::find(used_peak.begin(),used_peak.end(),double_hit_peaks[f])==used_peak.end()){
 std::cout<<"adding peak at bin #"<<double_hit_peaks[f]<<std::endl;
  FinalPeaks.push_back(double_hit_peaks[f]);
  MaxStartPoint.push_back(double_hit_peaks_start_point[f]);
  MaxEndPoint.push_back(double_hit_peaks_end_point[f]);
 }
}




 
 //..............END OF GROUPING 2-HIT-CLUSTERS
 
 



std::cout<<"---------------------*******-----------------------------------"<<std::endl;
std::cout<<"No of FinalPeaks after 2-hit evaluation = "<<FinalPeaks.size()<<std::endl;

for(unsigned int peak=0;peak<FinalPeaks.size();peak++){
std::cout<<"peak at bin # "<<FinalPeaks[peak]<<" ("<<-180+2*FinalPeaks[peak]<<" degrees). Its range is ["<<-180+2*MaxStartPoint[peak]<<", "<<-180+2*MaxEndPoint[peak]<<" ]"<<std::endl;
 }


std::cout<<"-----------------------*******---------------------------------"<<std::endl;

//.............. 2-hit peaks END............................



//Before I start actuall work on these 1-hit peaks I will get rid of all the ones that are too close to the already existing ranges:

std::vector<int> temp_one_hit_peaks;
temp_one_hit_peaks.clear();
int bad_one_hit_peak=0;

 for(unsigned int i=0;i<one_hit_peaks.size(); i++){
   for(unsigned int g=0; g<FinalPeaks.size(); g++){


   if(abs(one_hit_peaks[i]-MaxStartPoint[g])<6 || abs(one_hit_peaks[i]-MaxEndPoint[g])<6){
   
   bad_one_hit_peak=1;
   break;
   
   }
 }
 
 if (bad_one_hit_peak==0) temp_one_hit_peaks.push_back(one_hit_peaks[i]);
 bad_one_hit_peak=0;
}


one_hit_peaks.clear();
one_hit_peaks=temp_one_hit_peaks;


std::cout<<" ### NOW THE UPDATED NO OF one_hit_peaks= "<<one_hit_peaks.size()<<" they are at bins : ";

for(unsigned int u=0; u<one_hit_peaks.size();u++){

std::cout<<one_hit_peaks[u]<<std::endl;

}











//-----------WORK ON 1-HIT PEAKS----------------------------

std::vector<int> grouped_one_hit_peaks;
grouped_one_hit_peaks.clear();
int hand_made_peak=0;
int very_well_separated_two_hit_peak=0;
int two_hits_only=0;
int diff_between_one_hit_peaks=0;
int event_is_clean=0;

if(one_hit_peaks.size()<=4){event_is_clean=1;
std::cout<<"event_is_clean"<<std::endl;}

if(one_hit_peaks.size()>=4){ diff_between_one_hit_peaks=4;}

else{diff_between_one_hit_peaks=7;}

std::cout<<"one_hit_peaks size= "<<one_hit_peaks.size()<<std::endl;
if(one_hit_peaks.size()>=2){

  for(unsigned int i=0;i<one_hit_peaks.size(); i++){
    for(unsigned int j=0;j<one_hit_peaks.size(); j++){
    
      if(grouped_one_hit_peaks.size()==0){
      
      // loop thru existing ranges and make sure that your 1-hit peaks are not close to each other but also well separated from the existing ranges
       for(unsigned int r=0; r<FinalPeaks.size(); r++){ 
       std::cout<<"### FinalPeaks.size()= "<<FinalPeaks.size()<<std::endl;
       std::cout<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]<<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r]<<std::endl;
       if(i!=j && abs(one_hit_peaks[i]-one_hit_peaks[j])<=diff_between_one_hit_peaks ){ 
       
       //std::cout<<"abs(one_hit_peaks[i]-MaxStartPoint[r])= "<<abs(one_hit_peaks[i]-MaxStartPoint[r])<<" :: "<<one_hit_peaks[i]<<" - "<<MaxStartPoint[r]<<std::endl;
       
       
       //std::cout<<"MaxStartPoint["<<r<<"]= "<<MaxStartPoint[r]<<"MaxEndPoint["<<r<<"]= "<<MaxEndPoint[r]<<std::endl;
       
       //std::cout<<"one_hit_peaks[i]= "<<one_hit_peaks[i]<<" one_hit_peaks[j]= "<<one_hit_peaks[j]<<std::endl;
       
       
      //if you can't find this peak in the group, then add it:
         if(std::find(grouped_one_hit_peaks.begin(),grouped_one_hit_peaks.end(),one_hit_peaks[i])==grouped_one_hit_peaks.end()){
      grouped_one_hit_peaks.push_back(one_hit_peaks[i]);
      std::cout<<"ADDING PEAK #"<<one_hit_peaks[i]<<" TO THE GROUP"<<std::endl;
         }
        if(std::find(grouped_one_hit_peaks.begin(),grouped_one_hit_peaks.end(),one_hit_peaks[j])==grouped_one_hit_peaks.end()){
     grouped_one_hit_peaks.push_back(one_hit_peaks[j]);
     std::cout<<"ADDING PEAK #"<<one_hit_peaks[j]<<" TO THE GROUP"<<std::endl;
         }
   } //if difference between peak bin#s is small 
   }// loop thru already existing peaks 
} //if size=0 
    
     
    

     } //first loop thru each 1-hit peak
     
     
     
     
  } //second loop thru each 1-hit peak
  
  int failed=0;
   //Now go thru all the peaks again to find the third matching peak:
     if(grouped_one_hit_peaks.size()>0){
      for(unsigned int y=0;y<one_hit_peaks.size();y++){
       //now every added peak must be very close to EACH already added peak in the group
       for(unsigned int q=0; q<grouped_one_hit_peaks.size();q++){
       
        if(abs(one_hit_peaks[y]-grouped_one_hit_peaks[q])>5){
        failed=1;
        break;
        }
       
       } //for
       
       
       if(failed==0 && std::find(grouped_one_hit_peaks.begin(),grouped_one_hit_peaks.end(),one_hit_peaks[y])==grouped_one_hit_peaks.end() ){
       
       grouped_one_hit_peaks.push_back(one_hit_peaks[y]);
       
       //if(grouped_one_hit_peaks.size()==3) break;
       
       }
       failed=0;
       } //loop thru all 1-hit peaks
     }// if size >0
  
  
} //if we have more than 3 1-hit peaks


//so now we should have a group of 1-hit peaks that are close to each other. 

if(grouped_one_hit_peaks.size()>=3 || (grouped_one_hit_peaks.size()==2 && event_is_clean==1)){


 //form a peak by taking the average of them and picking the one that's closest to that value:
 int sum=0;
 int made_peak=200;
 
 
   for(unsigned int k=0; k<grouped_one_hit_peaks.size(); k++){
    
    sum+=grouped_one_hit_peaks[k];
  // std::cout<<"grouped_one_hit_peaks.size()= "<<grouped_one_hit_peaks.size()<<" sum= "<<sum<<std::endl;
   }
   
    made_peak=sum/grouped_one_hit_peaks.size();
    std::cout<<"---------MADE A HOME_MADE_PEAK FROM GROUPS OF 1-HIT CLUSTERS--------THIS PEAK IS AT BIN #"<<made_peak;
    FinalPeaks.push_back(made_peak);
    //Now work on the range. For now just start at the smallest peak bin-1 and end at the largest peak bin+1:
    
    std::sort(grouped_one_hit_peaks.begin(),grouped_one_hit_peaks.end());
    
    MaxStartPoint.push_back(grouped_one_hit_peaks[0]-1);
    MaxEndPoint.push_back(grouped_one_hit_peaks[grouped_one_hit_peaks.size()-1]+1);
    
    std::cout<<" its range is ["<<grouped_one_hit_peaks[0]-1<<", "<<grouped_one_hit_peaks[grouped_one_hit_peaks.size()-1]+1<<"]"<<std::endl;
    
    if(grouped_one_hit_peaks.size()>=3){hand_made_peak=1;}
    if(grouped_one_hit_peaks.size()==2){two_hits_only=1;}
    
}



//-----------END OF WORK ON 1-HIT PEAKS---------------------


std::cout<<"After BASIC EVALUATION we now have "<<FinalPeaks.size()<<" peaks"<<std::endl;
std::cout<<"FinalPeaks are at bin(s):  ";
for(unsigned int i=0; i<FinalPeaks.size();i++)
{
std::cout<<FinalPeaks[i]<<" which corresponds to angle=  "<<-180+2*FinalPeaks[i]<<std::endl;

}
//---------------------END BASIC EVALUATION-----------------------------

for(unsigned int peak=0; peak<MaxStartPoint.size(); peak++){


 //let's calculate how many bins away is the closest cluster to the one in question. You need to look to the left and to the right of each range to determine it. This is needed to pick the right MinHitsInRange value. Motivation: for very separated clusters in histos we want to be more lenient even though their peaks are not high. For very crowded environment want to be more strict.
 
  //loop thru all the other peaks to figure out the bin distance:
  for(unsigned int peak2=0; peak2<MaxStartPoint.size(); peak2++){
  
   if(peak!=peak2){
  diff_end_minus_start.push_back(MaxStartPoint[peak2]-MaxEndPoint[peak]);
  diff_start_minus_end.push_back(MaxStartPoint[peak]-MaxEndPoint[peak2]);
  
 
   }
  }
  
    for(unsigned int diff=0; diff<diff_end_minus_start.size();diff++){
    
    if(diff_end_minus_start[diff]>=0){
    positive_diff_end_minus_start.push_back(diff_end_minus_start[diff]); 
    
   
    }
    if(diff_start_minus_end[diff]>=0){
    positive_diff_start_minus_end.push_back(diff_start_minus_end[diff]); 
    
   
    }
    
    }
 //now take the minimum and this is your closest range in bin numbers:
 
 if(positive_diff_end_minus_start.size()>0){
 closest_range_right_side=*std::min_element(positive_diff_end_minus_start.begin(),positive_diff_end_minus_start.end());
 }
 else if(positive_diff_end_minus_start.size()==0){this_is_the_last_range=1;}
 
  if(positive_diff_start_minus_end.size()>0){
  closest_range_left_side=*std::min_element(positive_diff_start_minus_end.begin(),positive_diff_start_minus_end.end());
 }
 else if(positive_diff_start_minus_end.size()==0){this_is_the_first_range=1;}
 
 
 //if we only have 2 hits in the peak, it better be very well separated, otherwise not a valid peak!
 if((two_hits_only==1) && (closest_range_right_side>=16 || this_is_the_last_range==1 ) && (closest_range_left_side>=16 || this_is_the_first_range==1)){
 MinHitsInRange=2;
 very_well_separated_two_hit_peak=1;
 
 }
 
 
 
 
 
 
 
 
 
 
 
 
 //if range is well separated (or first or last):
 if((closest_range_right_side>=8 || this_is_the_last_range==1 ) && (closest_range_left_side>=8 || this_is_the_first_range==1)){
 MinHitsInRange=2;
 well_separated=1;
 std::cout<<"Peak at bin #"<<FinalPeaks[peak]<<" is ^^^ WELL SEPARATED ^^^"<<std::endl;
 
 if(this_is_the_last_range==1){std::cout<<" because it's the last range and ";}
 if(this_is_the_first_range==1){std::cout<<" because it's the first range and ";}
 if(closest_range_right_side>=8){std::cout<<"The closest range from right side is "<<closest_range_right_side<<" bins away"<<std::endl;}
 if(closest_range_left_side>=8){std::cout<<"The closest range from left side is "<<closest_range_left_side<<" bins away"<<std::endl;}
 
 }
 else{ MinHitsInRange=3; }
 
 
 //........................
 for(int bin=MaxStartPoint[peak]; bin<MaxEndPoint[peak];bin++){
 
 
 //std::cout<<"bin= "<<bin<<std::endl;
 if(plane==0){no_hits_in_range+=fh_theta_ind->GetBinContent(bin);
 //std::cout<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range<<std::endl;
 }
 if(plane==1){no_hits_in_range+=fh_theta_coll->GetBinContent(bin);
 //std::cout<<" plane= "<<plane<<" no_hits_in_range= "<<no_hits_in_range<<std::endl;
 }
 
 }
 
 
 
 
 std::cout<<"no_hits_in_range= "<<no_hits_in_range<<" for peak at bin # "<<FinalPeaks[peak]<<" ("<<-180+2*FinalPeaks[peak]<<" degrees). Its range is ["<<-180+2*MaxStartPoint[peak]<<", "<<-180+2*MaxEndPoint[peak]<<" ]"<<std::endl;
 
 if(plane==0){
 if((fh_theta_ind_Area->GetBinContent(FinalPeaks[peak])>0.4 && no_hits_in_range>=MinHitsInRange) || hand_made_peak==1 || very_well_separated_two_hit_peak==1){
 TempFinalPeaks.push_back(FinalPeaks[peak]);
 TempMaxStartPoint.push_back(MaxStartPoint[peak]);
 TempMaxEndPoint.push_back(MaxEndPoint[peak]);
 }
 else if(fh_theta_ind_Area->GetBinContent(FinalPeaks[peak])<=0.4 && no_hits_in_range>=MinHitsInRange && well_separated==1){
 TempFinalPeaks.push_back(FinalPeaks[peak]);
 TempMaxStartPoint.push_back(MaxStartPoint[peak]);
 TempMaxEndPoint.push_back(MaxEndPoint[peak]);
 
 }
 
 }//plane 0
 
 if(plane==1){
 if((fh_theta_coll_Area->GetBinContent(FinalPeaks[peak])>0.4 && no_hits_in_range>=MinHitsInRange) || hand_made_peak==1){
 TempFinalPeaks.push_back(FinalPeaks[peak]);
 TempMaxStartPoint.push_back(MaxStartPoint[peak]);
 TempMaxEndPoint.push_back(MaxEndPoint[peak]);
 }
 else if(fh_theta_coll_Area->GetBinContent(FinalPeaks[peak])<0.4 && no_hits_in_range>=MinHitsInRange && well_separated==1){
 TempFinalPeaks.push_back(FinalPeaks[peak]);
 TempMaxStartPoint.push_back(MaxStartPoint[peak]);
 TempMaxEndPoint.push_back(MaxEndPoint[peak]);
 
 }
 
 }//plane 1
 
 
 well_separated=0;
 no_hits_in_range=0;
 positive_diff_end_minus_start.clear();
 positive_diff_start_minus_end.clear();
 diff_end_minus_start.clear();
 diff_start_minus_end.clear();
 this_is_the_first_range=0;
 this_is_the_last_range=0;
 
 
 
} //for each peak


FinalPeaks=TempFinalPeaks;
MaxStartPoint=TempMaxStartPoint;
MaxEndPoint=TempMaxEndPoint;



std::cout<<" NO OF FINALPEAKS ***AFTER*** EVALUATION IS: "<<FinalPeaks.size()<<std::endl;

// If you just get 1 peak (which most likely corresponds to the muon track in case of CCQE) then let's check if there is a possibility of forming another cluster. Allow this formation if the event is CLEAN (small number of 1-hit peaks, and overall small no of all found peaks in the histos). Look for a peak that is VERY FAR AWAY from the already formed peak(s). But make sure that we get at least 2-hit in this newely formed cluster. Allow to run this test for both 1-peak-hits only since we already checked for 2-hit-peak possibilities above.

if(FinalPeaks.size()==1 && event_is_clean==1){
std::cout<<"$$$$  will try to come up with another peak since the event is clean"<<std::endl;
int separation=0;
std::vector<int> separated_one_hit_peaks, close_and_separated_one_h_pk;
separated_one_hit_peaks.clear();
close_and_separated_one_h_pk.clear();



//now take all the 1 hit peaks and see if you can find two which are very far from the already formed peak:

for(unsigned int onepk=0; onepk<one_hit_peaks.size();onepk++){
 for(unsigned int fpk=0;fpk<FinalPeaks.size();fpk++){
 //first pick the ones which are very far from the already existing peak
 if(abs(one_hit_peaks[onepk]-FinalPeaks[fpk])>17){
 
 separated_one_hit_peaks.push_back(one_hit_peaks[onepk]);
 std::cout<<"adding "<<one_hit_peaks[onepk]<<" to separated_one_hit_peaks"<<std::endl;
 
 }
 
 
 }




}


//Now check if you can find 2 1-hit-peaks which are also not so far from each other:

for(unsigned int i=0; i<separated_one_hit_peaks.size(); i++){
 for(unsigned int j=0; j<separated_one_hit_peaks.size(); j++){
 
 if(i!=j && abs(separated_one_hit_peaks[i]-separated_one_hit_peaks[j])<11){
 
 separation=abs(separated_one_hit_peaks[i]-separated_one_hit_peaks[j]);
 close_and_separated_one_h_pk.push_back(separated_one_hit_peaks[i]);
 close_and_separated_one_h_pk.push_back(separated_one_hit_peaks[j]);
 break;
 
 }
 
 
 }

if(close_and_separated_one_h_pk.size()==2) break;

}

//form a peak out of 2 1-hit-peaks:

if(close_and_separated_one_h_pk.size()==2){

int sum=0; 

for(unsigned int y=0; y<close_and_separated_one_h_pk.size(); y++){

sum+=close_and_separated_one_h_pk[y];
//and the range should be the start point of the first hit, the end point should be the end point of the second hit, so you could sort them. But, actually the range doesnt really matter here, just make it a few bins to the right and left!

}

int made_peak=sum/close_and_separated_one_h_pk.size();
std::cout<<"---------------------------------------------------------"<<std::endl;
std::cout<<" WARNING:: WILL MAKE AN EXCEPTION FOR THIS EVENT (event is clean, only 1 peak originally found) AND CREATE A NEW PEAK AT BIN #"<<made_peak<<std::endl;

std::cout<<"---------------------------------------------------------"<<std::endl;

FinalPeaks.push_back(made_peak);
MaxStartPoint.push_back(made_peak-0.5*separation);
MaxEndPoint.push_back(made_peak+0.5*separation);


}//if we have 2 1-hit-clusters



}


std::cout<<" FinalPeaks.size()="<<FinalPeaks.size()<<std::endl;

//-----------------------------------------------------------------
//One last check to Make sure that the peaks I have now are not too close to each other, otherwise throw the smallest which is very close. Right now they are ordered based on their area height, so we should be dropping them from the back, ie the smallest. (This piece of code might be move more up if needed, perhaps at the top of the basic evaluation section)

TempFinalPeaks.clear();
TempMaxStartPoint.clear();
TempMaxEndPoint.clear();


for( int pk=FinalPeaks.size()-1; pk>=0; pk--){
  for(unsigned int pk2=0; pk2<FinalPeaks.size();pk2++){

   if((unsigned int)(pk)!=pk2 && abs(FinalPeaks[pk]-FinalPeaks[pk2])<3){
  
    for(unsigned int i=0; i<FinalPeaks.size(); i++){
     if(i!=(unsigned int)(pk)){ 
     TempFinalPeaks.push_back(FinalPeaks[i]);
     TempMaxStartPoint.push_back(MaxStartPoint[i]);
     TempMaxEndPoint.push_back(MaxEndPoint[i]);
     }
    }
  
  FinalPeaks.clear();
  MaxStartPoint.clear();
  MaxEndPoint.clear();
  FinalPeaks=TempFinalPeaks;
  MaxStartPoint=TempMaxStartPoint;
  MaxEndPoint=TempMaxEndPoint;
  TempFinalPeaks.clear();
  TempMaxStartPoint.clear();
  TempMaxEndPoint.clear();
  break;
 
  }
 }
}

//..............................................................


std::cout<<" After making sure that each peak is more than 3 bins away from the previous one, FinalPeaks.size()="<<FinalPeaks.size()<<std::endl;

//-----------------------------------------------------------------

//Lastly, let's make sure that if there are ranges of peaks that lay right next to each other, their peak signal is actually different and not too small:
// I already set that the minimum peak must be greater than 0.4 so here the check is for small peaks of that order up to ~0.6

std::vector<int> bad_small_peak,marked;
bad_small_peak.clear();
marked.clear();

for(unsigned int peak=0; peak<MaxStartPoint.size(); peak++){
  for(unsigned int peak2=0; peak2<MaxStartPoint.size(); peak2++){

if((plane==0 && peak!=peak2 && abs(MaxStartPoint[peak]-MaxEndPoint[peak2])<=1 && (fh_theta_ind_Area->GetBinContent(FinalPeaks[peak])<0.6 || fh_theta_ind_Area->GetBinContent(FinalPeaks[peak2])<0.6)) || (plane==1 && peak!=peak2 && abs(MaxStartPoint[peak]-MaxEndPoint[peak2])<=1 && (fh_theta_coll_Area->GetBinContent(FinalPeaks[peak])<0.6 || fh_theta_coll_Area->GetBinContent(FinalPeaks[peak2])<0.6))) {

//get rid of one of them
 if(std::find(marked.begin(), marked.end(),peak)==marked.end() && std::find(marked.begin(), marked.end(),peak2)==marked.end()){
 bad_small_peak.push_back(FinalPeaks[peak]);
 marked.push_back(FinalPeaks[peak]);
 marked.push_back(FinalPeaks[peak2]);

 }//if not analyzed already

}




  }
}

//now copy the right peaks, if we found any small peak ranges right next to each other:
if(bad_small_peak.size()>0){
std::cout<<" ATTENTION: WILL NEED TO DELETE PEAKS AT THE FOLLOWING BIN #s, b/c its range is right next to some other range and the peak signal is < 0.6 "<<std::endl;

for(unsigned int bin=0; bin<bad_small_peak.size(); bin++){
std::cout<<bin<<std::endl;
}

TempFinalPeaks.clear();
TempMaxStartPoint.clear();
TempMaxEndPoint.clear();

for(unsigned int i=0; i<FinalPeaks.size(); i++){
 if(std::find(bad_small_peak.begin(),bad_small_peak.end(),FinalPeaks[i])==bad_small_peak.end()){
     TempFinalPeaks.push_back(FinalPeaks[i]);
     TempMaxStartPoint.push_back(MaxStartPoint[i]);
     TempMaxEndPoint.push_back(MaxEndPoint[i]);
 
 }

}

  bad_small_peak.clear();
  marked.clear();

  FinalPeaks.clear();
  MaxStartPoint.clear();
  MaxEndPoint.clear();
  FinalPeaks=TempFinalPeaks;
  MaxStartPoint=TempMaxStartPoint;
  MaxEndPoint=TempMaxEndPoint;
  TempFinalPeaks.clear();
  TempMaxStartPoint.clear();
  TempMaxEndPoint.clear();





}//if need to make corrections

//-----------------------------------------------------------------


need_to_reassign_hitsIDs=0;

std::cout<<"FORMING CLUSTERS NOW :) "<<std::endl;
std::cout<<"FinalPeaks are at bin(s):  ";
for(unsigned int i=0; i<FinalPeaks.size();i++)
{
std::cout<<FinalPeaks[i]<<" which corresponds to angle=  "<<-180+2*FinalPeaks[i]<<std::endl;

}


art::ServiceHandle<geo::Geometry> geom;
unsigned int channel=0, w=0;
unsigned int p=plane;
 unsigned int t=tpc;


//int no_noise_hits=0;
std::cout<<"In FindClusters(int plane), we should be producing "<<MaxStartPoint.size()<<" clusters"<<std::endl;
double a_polar, b_polar,theta_polar;

//HitsWithClusterID.clear();

std::vector<double> DiffAngles;
DiffAngles.clear();



for(unsigned int i = 0; i< allhits.size(); ++i){

 channel=allhits[i]->Wire()->RawDigit()->Channel();
 geom->ChannelToWire(channel,t,p,w);
 
 int diff_w= w - fwire_vertex_reco[plane];
      b_polar = diff_w*0.4; /**in cm*/
 //b_polar = (w - fwire_vertex[plane])* 0.4; /**in cm*/
 a_polar = (allhits[i]->PeakTime() - ftime_vertex_reco[plane])* ftimetick *fdriftvelocity; /** in cm*/
 theta_polar = fabs(asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2)))); /**in rad*/
 theta_polar = 180*theta_polar/fpi; /** in deg*/
 
 if(b_polar==0 && a_polar==0){
     theta_polar = 90;/** in deg*/
      }
     else if(b_polar==0 && a_polar<0){
     theta_polar = 180;/** in deg*/
      }
     else if(b_polar==0 && a_polar>0){
     theta_polar = 0;/** in deg*/
      }
     else if(b_polar>0 && a_polar==0){
     theta_polar = 90;/** in deg*/
      }
      else if(b_polar<0 && a_polar==0){
     theta_polar = -90;/** in deg*/
      }
      
 else if(b_polar>0 && a_polar>0){
     theta_polar = 90-theta_polar;/** in deg*/
       /** in deg*/
      }
      else if(b_polar>0 && a_polar<0){
       theta_polar = 90+theta_polar;/** in deg*/
      }
      else if(b_polar<0 && a_polar>0){
       theta_polar = -(90-theta_polar);/** in deg*/
      }
      else if(b_polar<0 && a_polar<0){
       theta_polar = -(90+theta_polar);/** in deg*/
      }
      
 for(unsigned int ClusterNo=0; ClusterNo<MaxStartPoint.size();ClusterNo++){

   if(theta_polar>=(-180+2*MaxStartPoint[ClusterNo]) && theta_polar<=(-180+2*MaxEndPoint[ClusterNo])){
     //want to start counting from 1, O is reserved for hits that will be marked as noise
     HitsWithClusterID.push_back(ClusterNo+1);
     break;}
   else if(ClusterNo==MaxStartPoint.size()-1){
     //decide where noise hits go
     
     //  if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
     //    std::cout<<"Noise hit at w= "<<w<<" t= "<<allhits[i]->PeakTime()<<" with theta_polar= "<<theta_polar;
     //   // std::cout<<"FinalPeaks.size()= "<<FinalPeaks.size()<<std::endl;
     //   }
     for(unsigned int peakNo=0;peakNo<FinalPeaks.size();peakNo++){
       DiffAngles.push_back(fabs(-180+2*FinalPeaks[peakNo]-theta_polar));
       // if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
       //    std::cout<<"diff for peak "<<peakNo<<" is "<<fabs(-180+2*FinalPeaks[peakNo]-theta_polar)<<std::endl;
       //    }
     }
     //now take minimum of DiffAngles and find at which position it is at, this position corresponds to clusterNo +1 , because we don't want to mark hits with zero 
     
     int position=std::distance(DiffAngles.begin(),std::min_element(DiffAngles.begin(),DiffAngles.end()));
   
     HitsWithClusterID.push_back(position+1);
     // if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
     //    std::cout<<"  This hit is closest to cluster # "<<position+1<<std::endl;
     //    }
     //no_noise_hits++;
     DiffAngles.clear();
   }


 } //loop over all ranges


 } //allhits


//std::cout<<"In FindClusters(int plane), marked "<<no_noise_hits<<" hits as NOISE"<<std::endl;
//std::cout<<"HitsWithClusterID contains the following clusterIDs:"<<std::endl;
 
//for(int i=0; i<HitsWithClusterID.size();i++){

//std::cout<<HitsWithClusterID[i]<<"  ";


//}

//std::cout<<std::endl;

//..............................................................................
//Now let's pick a minimum number of hits that we require in a cluster, MinHitsInCluster. If just formed cluster containes less than the desired number than it need to be assigned to the nearest cluster according to distance. This is done by removing the peak that corresponded to that cluster and rerunning assignement of hits again. We want to remove all the wrong peaks at once so reassignement is done only once. Put all the wrong peaks into a vector:

std::vector<unsigned int> WrongPeakNo;
WrongPeakNo.clear();
std::vector<int> WireNo;
int span=0;

int MinHitsInCluster=3; //later make it a parameter

for(unsigned int NClus=0; NClus<MaxStartPoint.size(); NClus++){
//search for clusters with too little hits (ie 1 or less than your desired parameter):

int NoHitsInCluster= std::count(HitsWithClusterID.begin(),HitsWithClusterID.end(),NClus+1);
std::cout<<"*** No of Hits for cluster # "<<NClus+1<<" = "<<NoHitsInCluster<<std::endl;

//........Check the span for small clusters................
  if(NoHitsInCluster<5 && NoHitsInCluster>0){
    
    WireNo.clear();
    for(unsigned int h=0; h<HitsWithClusterID.size();h++) {
    
        if(HitsWithClusterID[h]==NClus+1){
        
        WireNo.push_back(allhits[h]->Wire()->RawDigit()->Channel());
        }
    
    }

//now order WireNo and subtract first and last element to get the span in wire number:

std::sort(WireNo.begin(),WireNo.end());

span=WireNo[WireNo.size()-1]-WireNo[0];

  if(span>(NoHitsInCluster+0.66*NoHitsInCluster)){
    need_to_reassign_hitsIDs=1;
    WrongPeakNo.push_back(NClus);
  }


  }//if clusters containing less than 5 hits

//......end of checking span for small clusters...............
//
// let's give a chance for 2-hit-clusters IF it is very well separated from the rest of the peaks, say about 60 degrees which is 30 bins then do NOT reassign hits but leave this small 2-hit cluster as a valid one

 


  if(NoHitsInCluster<MinHitsInCluster && NoHitsInCluster!=2)
  {
  need_to_reassign_hitsIDs=1;


WrongPeakNo.push_back(NClus);
  }
  
  
  if(NoHitsInCluster==2)
  {
   for(unsigned int pk=0; pk<FinalPeaks.size();pk++){
    if(NClus!=pk && abs(FinalPeaks[NClus]-FinalPeaks[pk])<30){ go_ahead_at_reassign=1;
    break;}
   }
   
   if(go_ahead_at_reassign==1){WrongPeakNo.push_back(NClus);}
   if(go_ahead_at_reassign==0){
   
   std::cout<<"FOUND A 2-HIT-CLUSTER WHICH IS VERY WELL SEPARATED FROM THE REST, I will leave it (This is peak at bin #"<<FinalPeaks[NClus]<<std::endl;
   } 
   
  } //2-hit-clusters
  go_ahead_at_reassign=0;

}//loop thru all the clusters

if(need_to_reassign_hitsIDs==1){
//Now we need to change all the vectors containing info about peaks and subtract all the fake peaks:
HitsWithClusterID.clear();
std::vector<int> FinalPeaksTemporary;
std::vector<int> MaxStartPointTemporary;
std::vector<int> MaxEndPointTemporary;
FinalPeaksTemporary.clear();
MaxStartPointTemporary.clear();
MaxEndPointTemporary.clear();

std::vector<unsigned int>::iterator iter;
int initialNoPeaks=FinalPeaks.size();
for(int f=0; f<initialNoPeaks; f++){
iter=find(WrongPeakNo.begin(),WrongPeakNo.end(),f);
if(iter==WrongPeakNo.end()){
FinalPeaksTemporary.push_back(FinalPeaks[f]);
MaxStartPointTemporary.push_back(MaxStartPoint[f]);
MaxEndPointTemporary.push_back(MaxEndPoint[f]);

}

}
FinalPeaks.clear();
MaxStartPoint.clear();
MaxEndPoint.clear();

FinalPeaks=FinalPeaksTemporary;
MaxStartPoint=MaxStartPointTemporary;
MaxEndPoint=MaxEndPointTemporary;


}//if need to reassign
}













