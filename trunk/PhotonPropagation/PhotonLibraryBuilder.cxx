////////////////////////////////////////////////////////////////////////
/// \file    PhotonLibraryBuilder.cxx
/// \version $Id: CalHit.cxx,v 1.9 2011/01/24 23:18:19 p-novaart Exp $
/// \author  bpjones
////////////////////////////////////////////////////////////////////////
// PhotonLibraryBuilder.cxx  - Ben Jones, MIT 2010


#ifndef __CINT__

// LArSoft includes
#include "PhotonPropagation/PhotonLibraryBuilder.h"
#include "PhotonPropagation/PhotonLibraryService.h"
#include "Simulation/PhotonLibraryParameters.h"
#include "Simulation/PMTHit.h"

// FMWK includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Simulation/SimListUtils.h"

// ROOT includes
#include <TH1D.h>
#include <TF1.h>
#include <TTree.h>

// C++ language includes
#include <iostream>
#include <sstream>
#include <cstring>
#include <vector>

// Debug flag; only used during code development.
const bool debug = true;

namespace phot {
  

  PhotonLibraryBuilder::PhotonLibraryBuilder(fhicl::ParameterSet const& pset)
  {
    fLArG4InputModule = pset.get<std::string>("LArG4InputModule");
    fEvGenInputModule = pset.get<std::string>("EvGenInputModule");
  }
  

  void PhotonLibraryBuilder::beginJob()
  {
    art::ServiceHandle<PhotonLibraryService> pls;
    fLibrary = pls->GetLibrary();
  }

  
  PhotonLibraryBuilder::~PhotonLibraryBuilder() 
  {
  }
  
  void PhotonLibraryBuilder::analyze(art::Event const& evt)
  {
    //Get PMTHitCollection from Event
    art::ServiceHandle<sim::SimListUtils> slu;
    sim::PMTHitCollection TheHitCollection = slu->GetPMTHitCollection();

    art::Handle< std::vector<sim::PhotonLibraryParameters> > plHandle;
    evt.getByLabel(fEvGenInputModule,plHandle);
    art::Ptr<sim::PhotonLibraryParameters> ParamPtr;
    if(plHandle->size() == 1)    
      {
	art::Ptr<sim::PhotonLibraryParameters> p(plHandle, 0);
	ParamPtr=p;
      }
    else
      {
	std::cout<<"PhotonLibraryBuilder : No library parameters stored in event, cannot build library"<<std::endl;
	throw("");
      }
    
    if(ParamPtr->GetVoxelDef()!=fLibrary->GetVoxelDef())
      {
	std::cout<<"PhotonLibraryBuilder : Voxel definition of event generator and library do not match - cannot build library"<<std::endl;
	throw("");
      }
    
    int VoxID = ParamPtr->GetVoxelID();
    int Count = ParamPtr->GetPhotonCount();
    
    std::cout<<"PhotonLibraryBuilder : Adding entries to photon library for voxel " << VoxID<<std::endl;
    fLibrary->AddToHitCollectionForVoxel( VoxID,
					  TheHitCollection,
					  Count );
    
    
  }

  
  

  
}

#else
namespace sim
{
  class PhotonLibraryParameters;
}

#endif
