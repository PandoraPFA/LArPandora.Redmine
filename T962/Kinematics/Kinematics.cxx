////////////////////////////////////////////////////////////////////////
//
// VertexActivity class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to characterize the vertex activity associated with a neutrino event
////////////////////////////////////////////////////////////////////////


#include <iostream>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/OrphanHandle.h"


#include "T962/Kinematics/Kinematics.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"


#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "Simulation/LArVoxelCalculator.h"
#include "Simulation/LArVoxelData.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "T962/T962_Objects/ScanInfo.h"
#include "T962/T962_Objects/MINOS.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
kin::Kinematics::Kinematics(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel            (pset.get< std::string >("DBScanModuleLabel")),
  fLArG4ModuleLabel             (pset.get< std::string >("LArG4ModuleLabel")),
  fGenieGenModuleLabel             (pset.get< std::string >("GenieGenModuleLabel")),
  fScanModuleLabel              (pset.get< std::string > ("ScanModuleLabel")),
  fMINOSModuleLabel              (pset.get< std::string > ("MINOSModuleLabel"))
{
  //produces< std::vector<recob::Vertex> >();
}

//-----------------------------------------------------------------------------
kin::Kinematics::~Kinematics()
{
}

//-----------------------------------------------------------------------------


void kin::Kinematics::beginJob()
{
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;
  findcol    = tfs->make<TH2F>("IndCol",  "IndCol", 1000, -2, 500, 3000,-2,2048);
}

//-----------------------------------------------------------------------------
void kin::Kinematics::produce(art::Event& evt)
{
  double vertexcolwire=-1.;
  double vertexcoltime=-1.;
  double efficiency=1.;

  double elifetime_factor=0;
  double hitamplitude=0.;
  double vertex [3] = { 0, 0, 0 };
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larp;
  double electronlifetime=larp->ElectronLifetime();
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);

//hand scan info
   art::PtrVector<t962::ScanInfo> scanIn;
   scanIn.clear();
   art::Handle< std::vector<t962::ScanInfo> > scanHandle;
   evt.getByLabel(fScanModuleLabel,scanHandle);
 
  for(unsigned int i = 0; i < scanHandle->size(); ++i){     
      art::Ptr<t962::ScanInfo> scaninfo(scanHandle, i);
       scanIn.push_back(scaninfo);     
     }
 
     for(unsigned int i = 0; i < scanIn.size(); ++i){     
      vertexcolwire=scanIn[i]->Get_VertColWire();
      vertexcoltime=scanIn[i]->Get_VertColTime();   
     }
     
     std::cout<<"vertex: "<<vertexcolwire<<" "<<vertexcoltime<<std::endl;
     findcol->Fill(vertexcolwire,vertexcoltime);
     
     
 //minos info
   art::PtrVector<t962::MINOS> minosIn;
   minosIn.clear();
   art::Handle< std::vector<t962::MINOS> > minosHandle;
   evt.getByLabel(fMINOSModuleLabel,minosHandle);
 
  for(unsigned int i = 0; i < minosHandle->size(); ++i){     
      art::Ptr<t962::MINOS> minosinfo(minosHandle, i);
       minosIn.push_back(minosinfo);     
     }
 
     for(unsigned int i = 0; i < minosIn.size(); ++i){     
      std::cout<<"trkd cosy: "<<minosIn[i]->ftrkdcosy<<std::endl;
   
     }
//      
    
     
     
     
     
     
//take into account automated vertex finding as well (not yet implemented)  
//    art::PtrVector<recob::Vertex> vertIn;
//  
//    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
//      {
//        art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
//        vertIn.push_back(vertex);
//      }
//    
//     for(unsigned int i = 0; i < vertIn.size(); ++i){     
//       std::cout<<vertIn[i]->WireNum()<<" "<<vertIn[i]->DriftTime()<<std::endl;
//      }  

  filter::ChannelFilter chanFilt;  
  art::PtrVector<recob::Hit> cHits;
  art::PtrVector<recob::Hit> hit;
   
  art::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }




 }



