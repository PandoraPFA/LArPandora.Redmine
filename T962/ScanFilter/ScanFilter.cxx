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

#include "TMath.h"
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

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
    fMinShowers_req (pset.get< int >("MinShowers_req")),
    fMaxShowers_req (pset.get< int >("MaxShowers_req")),
    fMinTracks_req (pset.get< int >("MinTracks_req")),
    fMaxTracks_req (pset.get< int >("MaxTracks_req")),
    fFidVolume_cut (pset.get< float >("FidVolume_cut")),
    fCreateTTree   (pset.get< int >("CreateTTree")),
  fm_run(-1),
  fm_event(-1),
  fm_neutrino(-1),
  fm_maybeneutrino(-1),
  fm_tracks(-1),
  fm_showers(-1),
  fm_vertexx(-99.),
  fm_vertexy(-99.),
  fm_vertexz(-99.)
  {  
  }
  
  void ScanFilter::ScanFilter::beginJob()
{
  // get access to the TFile service
  if(fCreateTTree)
  {
  art::ServiceHandle<art::TFileService> tfs;
  ftree= tfs->make<TTree>("ScanTree","ScanTree");
  ftree->Branch("vertexx", &fm_vertexx , "vertexx/F");
  ftree->Branch("vertexy", &fm_vertexy , "vertexy/F");
  ftree->Branch("vertexz", &fm_vertexz , "vertexz/F");
  ftree->Branch("run", &fm_run , "run/I");
  ftree->Branch("event", &fm_event , "event/I");
  ftree->Branch("tracks", &fm_tracks , "tracks/I");
  ftree->Branch("showers", &fm_showers , "showers/I");
  ftree->Branch("neutrino", &fm_neutrino , "neutrino/I");
  ftree->Branch("maybeneutrino", &fm_maybeneutrino , "maybeneutrino/I");
  }
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
    float time=0;
    double y_vert=-99.;
    double z_vert=-99.;
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
 
 if(scanIn[i]->Get_IsnotNeutrino()!=1)
 {
    time=(scanIn[i]->Get_VertIndTime()+scanIn[i]->Get_VertColTime())/2.;
    
    //ArgoNeuT: cathode is at 60 ticks, induction plane is at 1630 ticks
    
  //   if(time>(1630-(fFidVolume_cut/.0297))||time<(60+(fFidVolume_cut/.0297)))
//     failFlag==1;
    
    int wire_I=scanIn[i]->Get_VertIndWire();
    
    int wire_C=scanIn[i]->Get_VertColWire();    
    wire_I=geom->PlaneWireToChannel(0,wire_I);
    wire_C=geom->PlaneWireToChannel(1,wire_C);
    geom->ChannelsIntersect(wire_I,wire_C,y_vert,z_vert);
    
    
//     double xyzstar[3];
//     double xyzend[3];
//     geom->WireEndPoints(0,120,xyzstar,xyzend);
//     
//     
//     std::cout<<xyzstar[0]<<" "<<xyzend[0]<<" "<<xyzstar[1]<<" "<<xyzend[1]<<std::endl;
  
     if(fCreateTTree)
     {
     fm_run=evt.id().run();
     fm_event=evt.id().event(); 
     fm_vertexx=(1630.-time)*.0297;
     fm_vertexy=y_vert;
     fm_vertexz=z_vert;
     fm_neutrino=scanIn[i]->Get_IsNeutrino();
     fm_maybeneutrino=scanIn[i]->Get_IsMaybeNeutrino();
     fm_tracks=scanIn[i]->Get_Track();
     fm_showers=scanIn[i]->Get_NumShower();     
     ftree->Fill();
     }
}
    if(    
    (scanIn[i]->Get_IsNeutrino()||(fNeutrino_req==1 && scanIn[i]->Get_IsMaybeNeutrino()==1)||(fNeutrino_req==0)) 
    && scanIn[i]->Get_NumShower()<=fMaxShowers_req
    && scanIn[i]->Get_NumShower()>=fMinShowers_req
    && scanIn[i]->Get_Track()<=fMaxTracks_req
    && scanIn[i]->Get_Track()>=fMinTracks_req
    && scanIn[i]->Get_Run()==run 
    && scanIn[i]->Get_Event()==event
 
    )
    
    failFlag=0;       
    }
 
    if(failFlag>0)
    return false;

    return true;
  }
	      
} //end namespace
