////////////////////////////////////////////////////////////////////////
//
// MergeData_Paddles class:
// For each ArgoNeuT event, find corresponding Paddles info. and 
// create a Paddles object to store the info. in the event record
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef MERGEDATAPADDLES_H
#define MERGEDATAPADDLES_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "T962/T962_Objects/Paddles.h"
#include "RawData/DAQHeader.h"
#include "RawData/BeamInfo.h"

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
#include <time.h>

///T962 analysis for merging different types of data
namespace merge {

  class MergeDataPaddles : public art::EDProducer {

  public:
          
    explicit MergeDataPaddles(fhicl::ParameterSet const& pset); 
    virtual ~MergeDataPaddles();

    void produce(art::Event& evt);
  
  private: 

    bool MergePaddles(t962::Paddles& paddles, art::Event& evt);   ///method to merge paddles data
    
    art::Handle<raw::DAQHeader> fdaq;
    std::string  fdaq_modulelabel;               ///< folder for input 
     
    
  }; // class MergeDataPaddles

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
      std::unique_ptr<t962::Paddles> Paddles(new t962::Paddles(paddles));
      evt.put(std::move(Paddles));
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

    std::ifstream myfile("lastevent.txt");
    std::string number;
    std::vector<std::string> stream_vec;
    if(myfile.is_open()){
      getline(myfile,number);
      
      std::istringstream stream;
      stream.clear();
      stream.str(number);
      stream_vec.clear();
      std::string times;
      while(std::getline(stream,times,' ')){
        if(!times.empty()) stream_vec.push_back(times);
      }
      
    }
    myfile.close();
    
    if(stream_vec.size()==0) std::cout << "Empty stream?" << std::endl;
    time_t lastdaqtime = (time_t) atoi(stream_vec[0].c_str());
    time_t lasttime = (time_t) atoi(stream_vec[1].c_str());
    
    std::cout << "daq time = " << spilltime << std::endl;
    std::cout << "lastdaqtime = " << lastdaqtime << std::endl;
    std::cout << "lasttime = " << lasttime << std::endl;
    
    bool foundpaddlesinfo = false;
    time_t time;
    std::vector<time_t> alltimes;
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
    int n = 0;
    int n_match = 0;
    int n_lastmatch = 0;
    
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
         alltimes.push_back(time);
         ++n;
         if(time==lasttime) n_lastmatch = n;
       }
       
       double diff = fabs(time - tms/1000.0);
       double diffc = (time - tms/1000.0);
       if(diff<10 && sv[0].compare(str0)==0)
         std::cout << "time = " << time << " tms = " << std::setprecision(13) << tms
                   << " diff = " << diffc << std::endl;
       
       if( diff<2.0 && time!=lasttime &&
           ( (n - n_lastmatch)==1
             || n_lastmatch==0 )){
          if(n_lastmatch==0 && mindiff!=1000.0 && sv[0].compare(str0)==0)
          {
             std::cout << "WARNING!  First event of run file could match to 2 Paddles events!" << std::endl;
             continue;
          }
         if(diff<mindiff){
           mindiff = diff;
           n_match = n;
           std::cout << "Setting mindiff...n = " << n << " n_lastmatch = " << n_lastmatch << std::endl;
         }
         foundpaddlesinfo = true;
       }
       if(diff<2.0 && time==lasttime){
          if(diff < mindiff && sv[0].compare(str0)==0) std::cout << "!!!!Event wants to match to last matched Paddles!" << std::endl;
       }
       if(sv[0].compare(str0)==0 && foundpaddlesinfo && diff==mindiff) paddles.SetTime(time);
       if(sv[0].compare(str1)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt1[i] = atoi(sv[i+1].c_str());
       if(sv[0].compare(str2)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt2[i] = atoi(sv[i+1].c_str());
       if(sv[0].compare(str3)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt3[i] = atoi(sv[i+1].c_str());
       if(sv[0].compare(str4)==0 && foundpaddlesinfo && diff==mindiff) for(int i = 0; i<4; ++i) pmt4[i] = atoi(sv[i+1].c_str());
   
       
       if(foundpaddlesinfo && diff>2.0){
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
      std::cout << "n_match = " << n_match << " n_lastmatch = " << n_lastmatch << std::endl;
      std::cout << paddles << std::endl;

      std::ofstream myfile;
      myfile.open("lastevent.txt");
      myfile << spilltime << " " << paddles.gettime() ;
      myfile.close();
      if((n_match - n_lastmatch)>1 && n_lastmatch!=0)
        std::cout << "Skipped Paddle event?" << std::endl;
    }
    else{
      std::cout << "No matching Paddle information found." << std::endl;
    }
   
  
    paddlesfile.close();
    
    return foundpaddlesinfo;
  
  }

  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(MergeDataPaddles);

}

#endif // MERGEDATAPADDLES_H
