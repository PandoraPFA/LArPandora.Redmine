////////////////////////////////////////////////////////////////////////
//
// DBSCANfinder.cxx
//
// kinga.partyka@yale.edu
//
//  This algorithm finds clusters of hits, they can be of arbitrary shape.You need to specify 2(3) parameters: 
// epsilon, epsilon2 and MinPoints as explained in the corresponding xml file.In my comments a 'point' reference 
// appears quite often. A 'point' is basically a simple hit which only contains wire and time information. This 
// algorithm is based on DBSCAN(Density Based Spatial Clustering of Applications with Noise): M. Ester, H.-P. Kriegel, 
// J. Sander, and X. Xu, A density-based algorithm for discovering clusters in large spatial databases with noise, 
// Second International Conference on Knowledge Discovery and Data Mining, pp. 226-231, AAAI Press. 1996. 
// ( Some of this code is from "Antonio Gulli's coding playground")  
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "ClusterFinder/DBScanService.h"
#include "ClusterFinder/DBcluster.h"
#include "Filters/ChannelFilter.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TGeoManager.h"
#include "TH1.h"

//-------------------------------------------------
cluster::DBcluster::DBcluster(fhicl::ParameterSet const& pset)
   
{  
  this->reconfigure(pset);
  produces<std::vector<recob::Cluster> >();  
}

//-------------------------------------------------
cluster::DBcluster::~DBcluster()
{
}

//-------------------------------------------------
void cluster::DBcluster::reconfigure(fhicl::ParameterSet const& p)
{
  fhitsModuleLabel=p.get< std::string >("HitsModuleLabel");
}

//-------------------------------------------------
void cluster::DBcluster::beginJob(){
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
  fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
  fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
    
}

//-----------------------------------------------------------------
void cluster::DBcluster::produce(art::Event& evt)
{
   
  //get a collection of clusters   
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    
  //std::cout << "event  : " << evt.Header().Event() << std::endl;
  art::ServiceHandle<geo::Geometry> geom;

  art::Handle< std::vector<recob::Hit> > hitcol;
  evt.getByLabel(fhitsModuleLabel,hitcol);
  
  
  ///loop over all hits in the event and look for clusters (for each plane)
  
  
  art::PtrVector<recob::Hit> allhits;
  art::PtrVector<recob::Hit> clusterHits;

  // get the DBScan service
  art::ServiceHandle<cluster::DBScanService> dbscan;
 
  // get the ChannelFilter
  filter::ChannelFilter chanFilt;
      
  unsigned int p(0),w(0), t(0), channel(0);
  for(unsigned int tpc = 0; tpc < geom->NTPC(); ++tpc){
    for(unsigned int plane=0; plane<geom->Nplanes(tpc); plane++){
      for(unsigned int i = 0; i< hitcol->size(); ++i){
      
	art::Ptr<recob::Hit> hit(hitcol, i);
      
	channel=hit->Wire()->RawDigit()->Channel();
	geom->ChannelToWire(channel,t,p,w);
    
	if(p == plane && t == tpc) allhits.push_back(hit);
	
      }  
      
      dbscan->InitScan(allhits, chanFilt.SetOfBadChannels());

      // std::cout<<"number of hits is: "<<hit.size()<<std::endl;
 
      //----------------------------------------------------------------
      for(unsigned int j = 0; j < dbscan->fps.size(); ++j){

	if(allhits.size() != dbscan->fps.size()) break;
   
	fhitwidth->Fill(dbscan->fps[j][2]);
	if(allhits[j]->Wire()->RawDigit()->Channel()<240){ fhitwidth_ind_test->Fill(dbscan->fps[j][2]);}
	if(allhits[j]->Wire()->RawDigit()->Channel()>240){ fhitwidth_coll_test->Fill(dbscan->fps[j][2]);}
	// std::cout<<"Point "<<j<<"= ("<<allhits[j]->Wire()->RawDigit()->Channel()
	//<<", "<<(allhits[j]->StartTime()+allhits[j]->EndTime())/2.<<", "<<dbscan->fps[j][2]<<")"<<std::endl;
	//  std::cout << dbscan->fps[j][i] << ' ';
	//  std::cout << std::endl;
      
      }
      //*******************************************************************

      dbscan->computeSimilarity();
      dbscan->computeSimilarity2();
      dbscan->computeWidthFactor();
      dbscan->run_cluster();

      //std::cout<<clusters;
      //std::cout<<"DBSCAN found "<<dbscan->fclusters.size()<<" cluster(s)."<<std::endl;


      for(unsigned int i=0; i<dbscan->fclusters.size();i++){
	//recob::Cluster* reco_cl= new recob::Cluster();
	for(unsigned int j=0;j<dbscan->fpointId_to_clusterId.size();j++){
	  
	  if(dbscan->fpointId_to_clusterId[j]==(i+1)){
	    
	    // reco_cl->Add(allhits[j]);
	    clusterHits.push_back(allhits[j]);
	  }
       
	}
         
	////////
	if (clusterHits.size()>0){

	  /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	  unsigned int sw = 0;
	  unsigned int ew = 0;
	  geom->ChannelToWire(clusterHits[0]->Wire()->RawDigit()->Channel(), t, p, sw);
	  geom->ChannelToWire(clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel(), t, p, ew);
	  	      
	  recob::Cluster cluster(clusterHits, 
				 sw*1., 0.,
				 clusterHits[0]->PeakTime(), clusterHits[0]->SigmaPeakTime(),
				 ew*1., 0.,
				 clusterHits[clusterHits.size()-1]->PeakTime(), clusterHits[clusterHits.size()-1]->SigmaPeakTime(),
				 -999., 0., 
				 -999., 0.,
				 i);

	  ccol->push_back(cluster);
	  //std::cout<<"no of hits for this cluster is "<<clusterHits.size()<<std::endl;
	  clusterHits.clear();
	  //////
	}//end if clusterHits has at least one hit
   
      }//end loop over fclusters

      allhits.clear();
    }//end loop over planes
  }//end loop over tpcs

  std::sort(ccol->begin(),ccol->end());//sort before Putting

  mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
  mf::LogVerbatim("Summary") << "DBcluster Summary:";
  for(unsigned int i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;

  evt.put(ccol);
  return;
}
