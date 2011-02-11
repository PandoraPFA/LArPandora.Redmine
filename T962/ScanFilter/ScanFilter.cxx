////////////////////////////////////////////////////////////////////////
//
// ScanFilter class
//
// joshua.spitz@yale.edu
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

//Framework Includes
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


//Larsoft Includes
#include "T962/ScanFilter/ScanFilter.h"
#include "RecoBase/recobase.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"// 

namespace filt{

  //-------------------------------------------------
  ScanFilter::ScanFilter(fhicl::ParameterSet const & pset) : 
    fScanModuleLabel(pset.get< std::string > ("ScanModuleLabel")),
    fNeutrino_req   (pset.get< int >("Neutrino_req")),
    fNumShowers_req (pset.get< int >("NumShowers_req")),
    fNumTracks_req  (pset.get< int >("NumTracks_req"))
  {   
  }

  //-------------------------------------------------
  ScanFilter::~ScanFilter()
  {
  }

  //-------------------------------------------------
  bool ScanFilter::filter(art::Event &evt)
  { 

    int failFlag = 1;
    int run = evt.id().run();
    int event = evt.id().event();
    
    art::PtrVector<t962::ScanInfo> scanIn;
    scanIn.clear();

    art::ServiceHandle<geo::Geometry> geom;
    art::Handle< std::vector<t962::ScanInfo> > scanHandle;
    evt.getByLabel(fScanModuleLabel,scanHandle);

    for(unsigned int i = 0; i < scanHandle->size(); ++i){     
     art::Ptr<t962::ScanInfo> scaninfo(scanHandle, i);
      scanIn.push_back(scaninfo);     
    }

    for(unsigned int i = 0; i < scanIn.size(); ++i){

    if(    
    (scanIn[i]->Get_IsNeutrino()||(fNeutrino_req==1 && scanIn[i]->Get_IsMaybeNeutrino()==1)||(fNeutrino_req==0 && scanIn[i]->Get_IsnotNeutrino()==1)) 
    && scanIn[i]->Get_NumShower()<=fNumShowers_req  
    && scanIn[i]->Get_Track()<=fNumTracks_req 
    && scanIn[i]->Get_Run()==run 
    && scanIn[i]->Get_Event()==event)
    failFlag=0;       
    }
 
    if(failFlag>0)
    return false;

    return true;
  }
	      
} //end namespace
