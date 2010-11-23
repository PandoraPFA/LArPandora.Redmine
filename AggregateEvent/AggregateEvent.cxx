////////////////////////////////////////////////////////////////////////
//
// Aggregate class
//
//
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

// Framework includes
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/View.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 

// LArSoft includes
#include "AggregateEvent/AggregateEvent.h"
#include "AggregateEvent/AggregateVertex.h"
#include "AggregateEvent/AggVertex.h"
#include "AggregateEvent/AggTrack.h"
#include "AggregateEvent/AggEvent.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
//#include "RecoBase/Calorimetry.h"
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/geo.h"

// ROOT includes
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"





//-------------------------------------------------
aggr::AggregateEvent::AggregateEvent(edm::ParameterSet const& pset) : 
  fCalorimetryModuleLabel(pset.getParameter< std::string >("CalorimetryModuleLabel"))
{

  produces< std::vector<aggr::AggEvent> >();

}

//-------------------------------------------------
aggr::AggregateEvent::~AggregateEvent()
{
  // Called, apparently, at end of (Sub?)Run, not Event. Porque?
}


//------------------------------------------------------------------------------------//
void aggr::AggregateEvent::produce(edm::Event& evt, edm::EventSetup const&)
{ 

  std::auto_ptr<std::vector<aggr::AggEvent> > ecol(new std::vector<aggr::AggEvent>);

  // get the geometry
  edm::Service<geo::Geometry> geom;

  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);
  for(unsigned int ii = 0; ii < hitListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Hit> hit(hitListHandle, ii);
      hitlist.push_back(hit); // class member
    }
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusterlist.push_back(cluster); // class member
    }
  edm::Handle< std::vector<recob::Cluster> > hclusterListHandle;
  evt.getByLabel(fHCModuleLabel,hclusterListHandle);
  for(unsigned int ii = 0; ii < hclusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> hcluster(hclusterListHandle, ii);
      hclusterlist.push_back(hcluster); // class member
    }
  edm::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      vertexlist.push_back(vertex); // class member
    }
  edm::Handle< std::vector<recob::Track> > trackListHandle;
  evt.getByLabel(fT3DModuleLabel,trackListHandle);
  for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Track> track(trackListHandle, ii);
      tracklist.push_back(track); // class member
    }
  edm::Handle< std::vector<recob::Shower> > showerListHandle;
  evt.getByLabel(fShowerModuleLabel,showerListHandle);
  for(unsigned int ii = 0; ii < showerListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Shower> shower(showerListHandle, ii);
      showerlist.push_back(shower); // class member
    }



  runNum = evt.run();
  evtNum = evt.id().event();



  std::cout << "aggr::AggregateEvent: Number of vertices is " << vertexlist.size() << std::endl;
  std::cout << "aggr::AggregateEvent: Number of clusters is " << clusterlist.size() << std::endl;
  std::cout << "aggr::AggregateEvent: Number of Houghclusters is " << hclusterlist.size() << std::endl;
  std::cout << "aggr::AggregateEvent: Number of showers is " << showerlist.size() << std::endl;
  std::cout << "aggr::AggregateEvent: Number of 3D tracks is " << tracklist.size() << std::endl;




}





