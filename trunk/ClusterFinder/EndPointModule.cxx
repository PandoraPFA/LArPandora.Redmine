////////////////////////////////////////////////////////////////////////
//
// EndPointModule class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find (weak) vertices from hits after deconvolution and hit finding. 
//  A weak vertex is a vertex that has been found using a dedicated vertex finding algorithm only. A 
//  strong vertex is a vertex that has been found using a dedicated vertex finding algorithm and matched 
//  to a crossing of two or more HoughLineFinder lines. The VertexMatch module finds strong vertices.
////////////////////////////////////////////////////////////////////////
/// The algorithm is based on:
///C. Harris and M. Stephens (1988). "A combined corner and edge detector". Proceedings of the 4th Alvey 
///Vision Conference. pp. 147-151.
///B. Morgan (2010). "Interest Point Detection for Reconstruction in High Granularity Tracking Detectors". 
///arXiv:1006.3012v1 [physics.ins-det]
//Thanks to B. Morgan of U. of Warwick for comments and suggestions

#include <iostream>

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

#include "ClusterFinder/EndPointService.h"
#include "ClusterFinder/EndPointModule.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"


#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
cluster::EndPointModule::EndPointModule(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel  (pset.get< std::string >("DBScanModuleLabel"))
  
{
  produces< std::vector<recob::EndPoint2D> >();
}

//-----------------------------------------------------------------------------
cluster::EndPointModule::~EndPointModule()
{
}

//-----------------------------------------------------------------------------

void cluster::EndPointModule::produce(art::Event& evt)
{

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  //Point to a collection of vertices to output.
  
  //.......................................
  art::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
 
  art::ServiceHandle<cluster::EndPointService> vs;
  
  // make a std::vector<recob::Cluster> for the output of the 
  // Hough Transform
  std::vector<recob::EndPoint2D> vtxOut;
  
  size_t numvtx = vs->EndPoint(clusIn, vtxOut);

  LOG_DEBUG("Vertex") << "found " << numvtx << "vertices with VertexService";

  //Point to a collection of vertices to output.
  std::auto_ptr<std::vector<recob::EndPoint2D> > vtxcol(new std::vector<recob::EndPoint2D>(vtxOut));

  
  evt.put(vtxcol);   
}

