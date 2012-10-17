////////////////////////////////////////////////////////////////////////
//
// MergeData_Beam class:
// For each ArgoNeuT event, find corresponding NuMI spill info. and 
// create a BeamInfo object to store the info. in the event record
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGEDATABEAM_H
#define MERGEDATABEAM_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// LArSoft include
#include "RawData/BeamInfo.h"
#include "RawData/DAQHeader.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <dirent.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>



///T962 analysis for merging different types of data
namespace merge {

  class MergeDataBeam : public art::EDProducer {

  public:
          
    explicit MergeDataBeam(fhicl::ParameterSet const& pset); 
    virtual ~MergeDataBeam();

    void produce(art::Event& evt);
  
  private: 

    bool MergeBeam(raw::BeamInfo& beam);                 ///method to merge beam data
    
    art::Handle<raw::DAQHeader> fdaq;
    std::string  fdaq_modulelabel;               ///< folder for input 
     
    
  }; // class MergeDataBeam

  //-------------------------------------------------
  MergeDataBeam::MergeDataBeam(fhicl::ParameterSet const& pset) : 
    fdaq_modulelabel(pset.get< std::string >("daq"))
  {
    produces<raw::BeamInfo> ();
  }

  //-------------------------------------------------
  MergeDataBeam::~MergeDataBeam()
  {
    
  }

  //-------------------------------------------------
  void MergeDataBeam::produce(art::Event& evt)
  {
    evt.getByLabel(fdaq_modulelabel,fdaq);
       
    raw::BeamInfo beam;

    if(MergeBeam(beam)){
      std::unique_ptr<raw::BeamInfo> Beam(new raw::BeamInfo(beam));
      evt.put(std::move(Beam));
    }
  
    return;
  }
  
  //-------------------------------------------------
  bool MergeDataBeam::MergeBeam(raw::BeamInfo& beam)
  {
    time_t spilltime = fdaq->GetTimeStamp();//time info. from DAQ480 software
    spilltime = spilltime >> 32;
    tm *timeinfo = localtime(&spilltime);//use this to get day/month/year;we should worry about users in different time zones at some point
   
    char beamfilename[20];
    //sprintf(beamfilename,"/argoneut/data/rundata/matched/matched_%02d_%i_%d",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
    sprintf(beamfilename,"./matchfiles/matched_%02d_%i_%d",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
    std::ifstream beamfile(beamfilename);

    if(!beamfile.is_open()){
       LOG_DEBUG ("badfile") << "MergeBeam could not open file named " << beamfilename ;
       return false;
    }
  
    long long int tms;
    double tor101;
    double tortgt;
    double trtgtd;
    std::string first,k;
    int event,run;
  
    bool foundbeaminfo = false;
  
    while(getline(beamfile,k)){
    
       std::istringstream ins;
       ins.clear();
       ins.str(k);
       ins>>first>>run>>event>>tms>> tor101 >> tortgt >> trtgtd;
//        std::cout<<std::setprecision(9)<<first<<"  "
// 		<<run<<"  "<<event<<"  "<<tms<<" "
// 		<<tor101<<"  "<< tortgt <<"  "<< trtgtd <<std::endl;
   
       if((fdaq->GetRun()==run)&&(fdaq->GetEvent()==event ) ){ 
          foundbeaminfo=true;
          beam.SetTOR101(tor101);
          beam.SetTORTGT(tortgt);
          beam.SetTRTGTD(trtgtd);
          beam.SetT_MS(tms);
          break;
       }
    }
  
    beamfile.close();
    
    return foundbeaminfo;
  
  }

  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(MergeDataBeam);

}

#endif // MERGEDATABEAM_H
