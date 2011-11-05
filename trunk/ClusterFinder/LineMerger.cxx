////////////////////////////////////////////////////////////////////////
//
// LineMerger class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// biagio.rossi@lhep.unibe.ch
// msoderbe@syr.edu
// joshua.spitz@yale.edu
//
// This algorithm is designed to merge 2D lines with similar slope and endpoints 
//  
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

//LArSoft includes:
#include "ClusterFinder/LineMerger.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>


namespace cluster{

  //-------------------------------------------------
  LineMerger::LineMerger(fhicl::ParameterSet const& pset) : 
    fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel")),
    fSlope             (pset.get<double     >("Slope")),
    fEndpointWindow         (pset.get<double     >("EndpointWindow"))
  {
    produces< std::vector<recob::Cluster> >();
  }

  //-------------------------------------------------
  LineMerger::~LineMerger()
  {
  }

  //-------------------------------------------------
  void LineMerger::beginJob()
  {
    //this doesn't do anything now, but it might someday
  }
    
  //------------------------------------------------------------------------------------//
  void LineMerger::produce(art::Event& evt)
  { 
    // Get a Handle for the input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel,clusterVecHandle);

    art::ServiceHandle<geo::Geometry> geo;
    int nplanes = geo->Nplanes();

    art::PtrVector<recob::Cluster> Cls[nplanes];//one PtrVector for each plane in the geometry
    std::vector<int> Cls_matches[nplanes];//vector with indicators for whether a cluster has been merged already

    // loop over the input Clusters
    for(unsigned int i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a art::Ptr to each Cluster
      art::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      
      switch(cl->View()){
      case geo::kU :
	Cls[0].push_back(cl);
	Cls_matches[0].push_back(0);
	break;
      case geo::kV :
	Cls[1].push_back(cl);
	Cls_matches[1].push_back(0);
	break;
      case geo::kW :
	Cls[2].push_back(cl);
	Cls_matches[2].push_back(0);
	break;
      default :
	break;
      }
    }

     //////////////////////////////////////////////////////
    // Make a std::auto_ptr<> for the thing you want to put into the event
    // because that handles the memory management for you
    //////////////////////////////////////////////////////
    std::auto_ptr<std::vector<recob::Cluster> > SuperClusters(new std::vector<recob::Cluster>);

    for(int i = 0; i<nplanes; ++i){
      int clustersfound = 0;//how many merged clusters found in each plane
      int clsnum1 = 0;
      for(art::PtrVector<recob::Cluster>::const_iterator clusIter1 = Cls[i].begin(); clusIter1!=Cls[i].end();++clusIter1){
	art::Ptr<recob::Cluster> cl1 = *clusIter1;
	if(Cls_matches[i][clsnum1]==1){
	  clsnum1++;
	  continue;
	}
	SuperClusters->push_back(*cl1); 
	Cls_matches[i][clsnum1]=1; 
	//SuperClusters->back().SetID(clustersfound);//IDs are sequential by plane, starting from 0
	++clustersfound;
	// recob::Cluster SCl= SuperClusters->back(); //spitz commented this out and moved it inside the for loop.
	
	int clsnum2 = 0;
	for(art::PtrVector<recob::Cluster>::const_iterator clusIter2 = Cls[i].begin(); clusIter2!=Cls[i].end();++clusIter2){
	  art::Ptr<recob::Cluster> cl2 = *clusIter2;
	  recob::Cluster SCl= SuperClusters->back(); //spitz moved this inside the for loop.
	  if(Cls_matches[i][clsnum2]==1){
	    clsnum2++;
	    continue;
	  }

	  //check that the slopes are the same
	  //added 13.5 ticks/wirelength in ArgoNeuT. need to make this detector agnostic--spitz
	  //would be nice to have a LArProperties function that returns ticks/wire.
	  bool sameSlope = SlopeCompatibility(SCl.dTdW()*(1./13.5),cl2->dTdW()*(1./13.5));  

	  //check that the endpoints fall within a circular window of each other //spitz did this in place of intercept matching
	  bool sameEndpoint = EndpointCompatibility(SCl.StartPos(),SCl.EndPos(),cl2->StartPos(),cl2->EndPos());


	  if(sameSlope && sameEndpoint){
	    SuperClusters->back() = SuperClusters->back() + *cl2;
	    Cls_matches[i][clsnum2]=1;       
	  }
	  clsnum2++;
	}
	clsnum1++;
      }
    }

    std::sort(SuperClusters->begin(),SuperClusters->end());//sort before Putting

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "LineMerger Summary:";
    for(unsigned int i = 0; i<SuperClusters->size(); ++i) mf::LogVerbatim("Summary") << SuperClusters->at(i) ;

    evt.put(SuperClusters);
     
    return;

  }

  //------------------------------------------------------------------------------------//
  //checks the difference between angles of the two lines
  bool LineMerger::SlopeCompatibility(double slope1, double slope2)
  { 
    double sl1=atan(slope1);
    double sl2=atan(slope2);
    bool comp = fabs(sl1-sl2)<fSlope ? true : false;//the units of fSlope are radians

    return comp;
  }
  //------------------------------------------------------------------------------------//
  bool LineMerger::EndpointCompatibility(std::vector<double> sclstart, std::vector<double> sclend,std::vector<double> cl2start, std::vector<double> cl2end)
  { 
   	  double sclstartwire=sclstart[0];
	  double sclstarttime=sclstart[1];
	  double sclendwire=sclend[0];
	  double sclendtime=sclend[1];
	  
	  double cl2startwire=cl2start[0];
	  double cl2starttime=cl2start[1];
	  double cl2endwire=cl2end[0];
	  double cl2endtime=cl2end[1];

	 //13.5 ticks/wire. need to make this detector agnostic--spitz
     double distance=sqrt((pow(sclendwire-cl2startwire,2)*13.5)+pow(sclendtime-cl2starttime,2));
     //not sure if this line is necessary--spitz
     double distance2=sqrt((pow(sclstartwire-cl2endwire,2)*13.5)+pow(sclstarttime-cl2endtime,2));

    bool comp = (distance<fEndpointWindow||distance2<fEndpointWindow) ? true : false;
    return comp;
  }

}//namespace cluster{
