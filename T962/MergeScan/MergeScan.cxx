////////////////////////////////////////////////////////////////////////
//
//   MergeScan Class
//
//    joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <math.h>
#include "T962/MergeScan/MergeScan.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include<dirent.h>
#include <time.h>
#include <stdio.h>
#include<cstdlib>

namespace t962{

//-------------------------------------------------
MergeScan::MergeScan(fhicl::ParameterSet const& pset) : 
  
  daq_modulelabel     (pset.get< std::string >("daq")),
  scanners            (pset.get< std::string >("scanners")),  
  foundscaninfo(false)
{
produces< std::vector<t962::ScanInfo> >();
}
  
void MergeScan::beginJob()
{

  // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
}

//-------------------------------------------------
MergeScan::~MergeScan()
{
 
}

//-------------------------------------------------

void MergeScan::produce(art::Event& evt)
{
  std::auto_ptr<std::vector<t962::ScanInfo> > Scan_coll(new std::vector<t962::ScanInfo> );
  //art::Handle< std::vector<raw::DAQHeader> > daqHandle;

  evt.getByLabel(daq_modulelabel,fdaq);
  //art::Ptr<raw::DAQHeader> daq = art::Ptr<raw::DAQHeader>(daqHandle, daqHandle->size()-1);
  t962::ScanInfo scan;

  time_t spilltime = fdaq->GetTimeStamp();//time info. from DAQ480 software
  tm *timeinfo = localtime(&spilltime);

  char scanfilename[100];
 sprintf(scanfilename,"/argoneut/data/simplescan_data/simplescan_text/scan_%s.00%d.txt",scanners.c_str(),fdaq->GetRun());

 std::ifstream scanfile(scanfilename);

  if(!scanfile.is_open()){
    std::cerr << "MergeScan:  Could not open file named " << scanfilename << std::endl;
    return;
  }

  int isneutrino,isnotneutrino,ismaybeneutrino,track,vertindtime,vertcoltime,vertindwire,vertcolwire,numshower;
  std::string first,k;
  int event,run;

  while(getline(scanfile,k)){

    std::istringstream ins3;
    ins3.clear();
    ins3.str(k);
        ins3>>run>>event>>isnotneutrino>>ismaybeneutrino>>isneutrino>>track>>vertindtime>>vertcoltime>>vertindwire>>vertcolwire>>numshower;
        
    if((fdaq->GetRun()==run)&&(fdaq->GetEvent()==event ) )
    foundscaninfo=true; 
    else
    foundscaninfo=false;
        
    if(foundscaninfo)
    {  
      scan.SetRun(run);
      scan.SetEvent(event);
      scan.SetIsNeutrino(isneutrino);
      scan.SetIsnotNeutrino(isnotneutrino);
      scan.SetIsMaybeNeutrino(ismaybeneutrino);
      scan.SetTrack(track);     
      scan.SetVertIndTime(vertindtime);
      scan.SetVertColTime(vertcoltime);
      scan.SetVertIndWire(vertindwire);
      scan.SetVertColWire(vertcolwire);
      scan.SetNumShower(numshower);
      // scan.SetScanner(scanners);
      break;
    }
  }
  scanfile.close();
  
  if(foundscaninfo)
  Scan_coll->push_back(scan);
  
  evt.put(Scan_coll);
  return;
  
}

}
