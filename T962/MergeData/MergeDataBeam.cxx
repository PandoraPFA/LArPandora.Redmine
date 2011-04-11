////////////////////////////////////////////////////////////////////////
//
//   MergeData_Beam Class
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <dirent.h>
#include <time.h>
#include <stdio.h>
#include <cstdlib>

#include "T962/MergeData/MergeDataBeam.h"

namespace merge{

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
      std::auto_ptr<raw::BeamInfo> Beam(new raw::BeamInfo(beam));
      evt.put(Beam);
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
    sprintf(beamfilename,"/argoneut/data/rundata/matched/matched_%02d_%i_%d",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
    std::ifstream beamfile(beamfilename);

    if(!beamfile.is_open()){
       mf::LogDebug ("badfile") << "MergeBeam could not open file named " << beamfilename ;
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

}//namespace merge
