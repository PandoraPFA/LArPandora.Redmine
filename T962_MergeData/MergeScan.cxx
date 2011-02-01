////////////////////////////////////////////////////////////////////////
//
//   MergeScan Class
//
//    joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <math.h>
#include "T962_MergeData/MergeScan.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include<dirent.h>
#include <time.h>
#include <stdio.h>
#include<cstdlib>

namespace merge{

//-------------------------------------------------
MergeScan::MergeScan(fhicl::ParameterSet const& pset) : 
  
  daq_modulelabel     (pset.get< std::string >("daq")),
  scanners            (pset.get< std::vector<std::string> >("scanners")),  
  foundscaninfo(false)
{
produces< std::vector<merge::ScanInfo> >();
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
  std::auto_ptr<std::vector<merge::ScanInfo> > Scan_coll(new std::vector<merge::ScanInfo> );
  art::Handle< std::vector<raw::DAQHeader> > daqHandle;
  evt.getByLabel(daq_modulelabel,daqHandle);
  art::Ptr<raw::DAQHeader> daq = art::Ptr<raw::DAQHeader>(daqHandle, daqHandle->size()-1);
  merge::ScanInfo scan;

  time_t spilltime = daq->GetTimeStamp();//time info. from DAQ480 software
  tm *timeinfo = localtime(&spilltime);
  for(int scanner=0;scanner<scanners.size();scanner++)
  {
  char scanfilename[100];
 sprintf(scanfilename,"/argoneut/data/simplescan_data/simplescan_text/scan_%s.00%d.txt",scanners[scanner].c_str(),daq->GetRun());

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
        
    if((daq->GetRun()==run)&&(daq->GetEvent()==event ) )
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
      scan.SetScanner(scanner);
      break;
    }
  }
  scanfile.close();
  
  if(foundscaninfo)
  Scan_coll->push_back(scan);
  }
  evt.put(Scan_coll);
  return;
  
}

}
