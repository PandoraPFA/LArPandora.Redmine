////////////////////////////////////////////////////////////////////////
//
//   MergeData_Paddles Class
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
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
#include <time.h>

#include "T962/MergeData/MergeDataPaddles.h"
#include "RawData/BeamInfo.h"

namespace merge{

  //-------------------------------------------------
  MergeDataPaddles::MergeDataPaddles(fhicl::ParameterSet const& pset) : 
    fdaq_modulelabel(pset.get< std::string >("daq"))
  {
    produces<t962::Paddles> ();
  }

  //-------------------------------------------------
  MergeDataPaddles::~MergeDataPaddles()
  {
    
  }

  //-------------------------------------------------
  void MergeDataPaddles::produce(art::Event& evt)
  {
    evt.getByLabel(fdaq_modulelabel,fdaq);
       
    t962::Paddles paddles;

    if(MergePaddles(paddles,evt)){
      std::auto_ptr<t962::Paddles> Paddles(new t962::Paddles(paddles));
      evt.put(Paddles);
    }
  
    return;
  }
  
  //-------------------------------------------------
  bool MergeDataPaddles::MergePaddles(t962::Paddles& paddles, art::Event& evt)
  {
   
    time_t spilltime = fdaq->GetTimeStamp();//time info. from DAQ480 software
    spilltime = spilltime >> 32;
    tm *timeinfo = localtime(&spilltime);//use this to get day/month/year;we should worry about users in different time zones at some point
    
    char paddlesfilename[100];
    sprintf(paddlesfilename,"/argoneut/data/paddles/pmt_%02d_%i_%d.txt",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
    
    //printf("Opening /argoneut/data/paddles/pmt_%02d_%i_%d.txt",timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year+1900);
    std::ifstream paddlesfile(paddlesfilename);
    
    if(!paddlesfile.is_open()){
       LOG_DEBUG ("badfile") << "MergePaddles could not open file named " << paddlesfilename ;
       return false;
    }
    
    bool foundpaddlesinfo = false;
    time_t time;
    int pmt1[4] = {0};
    int pmt2[4] = {0};
    int pmt3[4] = {0};
    int pmt4[4] = {0};

    std::string str0 ("time");
    std::string str1 ("pmt1");
    std::string str2 ("pmt2");
    std::string str3 ("pmt3");
    std::string str4 ("pmt4");
 
    std::string k;

    art::Handle< raw::BeamInfo > beam;
    evt.getByLabel("beam",beam);
    double tms = (double)beam->get_t_ms();
    double mindiff = 1000.0;
    
    while(getline(paddlesfile,k)){
    
       std::istringstream ins;
       ins.clear();
       ins.str(k);
       std::vector<std::string> sv;
       sv.clear();
       std::string word;
       while(std::getline(ins, word,' ')){
         if(!word.empty()) sv.push_back(word);
       }
       if(sv.size()==0) continue;
       if(sv[0].compare(str0)==0){
         time = atoi(sv[1].c_str());
       }
       
       double diff = fabs(time - tms/1000.0);
     
       double diffc = (time - tms/1000.0);
     //   if(diff<10 && sv[0].compare(str0)==0)
//          std::cout << "time = " << time << " tms = " << std::setprecision(13) << tms
//                    << " diff = " << diffc << std::endl;
 
       if(diff<1.5 && diffc<0.0) {
         if(diff<mindiff) mindiff = diff;
         foundpaddlesinfo = true;
       }
       if(sv[0].compare(str0)==0 && foundpaddlesinfo && diff==mindiff) paddles.SetTime(time);
       if(sv[0].compare(str1)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt1[i] = atoi(sv[i+1].c_str());
       if(sv[0].compare(str2)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt2[i] = atoi(sv[i+1].c_str());
       if(sv[0].compare(str3)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt3[i] = atoi(sv[i+1].c_str());
       if(sv[0].compare(str4)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt4[i] = atoi(sv[i+1].c_str());
   
       
       if(foundpaddlesinfo && diff>1.5){
         paddles.SetPMT(0,pmt1);
         paddles.SetPMT(1,pmt2);
         paddles.SetPMT(2,pmt3);
         paddles.SetPMT(3,pmt4);
         break;
       }
    }

    
    if(foundpaddlesinfo){
      std::cout << "Matching Paddle information found!" << std::endl;
      std::cout << "beam tms = " << std::setprecision(13) << tms << std::endl;  
      std::cout << "daq time = " << spilltime << std::endl;
      std::cout << paddles << std::endl;
    }
    else{
      std::cout << "No matching Paddle information found." << std::endl;
    }
   
  
    paddlesfile.close();
    
    return foundpaddlesinfo;
  
  }

}//namespace merge
