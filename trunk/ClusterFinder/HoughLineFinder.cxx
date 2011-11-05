////////////////////////////////////////////////////////////////////////
//
// \file HoughLineFinder.cxx
//
// \author joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN 
//  after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine 
//  Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

// ROOT includes
#include <TCanvas.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

// ART includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes 
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "ClusterFinder/HoughLineService.h"
#include "ClusterFinder/HoughLineFinder.h"


cluster::HoughLineFinder::HoughLineFinder(fhicl::ParameterSet const& pset) : 
  fDBScanModuleLabel       (pset.get< std::string >("DBScanModuleLabel"))
  
{
  produces< std::vector<recob::Cluster> >();
}

cluster::HoughLineFinder::~HoughLineFinder()
{
}



void cluster::HoughLineFinder::produce(art::Event& evt)
{

  //////////////////////////////////////////////////////
  // here is how to get a collection of objects out of the file
  // and connect it to a art::Handle
  //////////////////////////////////////////////////////
  // Read in the clusterList object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);

  art::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
  
  art::ServiceHandle<cluster::HoughLineService> hls;
  
  // make a std::vector<recob::Cluster> for the output of the 
  // Hough Transform
  std::vector<recob::Cluster> clusOut;
  
  size_t numclus = hls->Transform(clusIn, clusOut);

  LOG_DEBUG("HoughLineClusters") << "found " << numclus << "clusters with HoughLineService";

  //Point to a collection of clusters to output.
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(clusOut));

  std::sort(ccol->begin(),ccol->end());//sort before Putting

  mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
  mf::LogVerbatim("Summary") << "HoughLineFinder Summary:";
  for(size_t i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;

  evt.put(ccol);
 
    

 
}



