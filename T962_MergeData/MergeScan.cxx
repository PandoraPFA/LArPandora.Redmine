////////////////////////////////////////////////////////////////////////
//
//   MergeScan Class
//
//    joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
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
MergeScan::MergeScan(edm::ParameterSet const& pset) : 
  
  daq_modulelabel(pset.getParameter< std::string >("daq")),  
  foundscaninfo(false)
{
produces< std::vector<merge::ScanInfo> >();
}
  
void MergeScan::beginJob(edm::EventSetup const&)
{

  // get access to the TFile service
    edm::Service<edm::TFileService> tfs;
}

//-------------------------------------------------
MergeScan::~MergeScan()
{
 
}

//-------------------------------------------------

void MergeScan::produce(edm::Event& evt, edm::EventSetup const&)
{
  std::auto_ptr<std::vector<merge::ScanInfo> > Scan_coll(new std::vector<merge::ScanInfo> );
  edm::Handle< std::vector<raw::DAQHeader> > daqHandle;
  evt.getByLabel(daq_modulelabel,daqHandle);
  edm::Ptr<raw::DAQHeader> daq = edm::Ptr<raw::DAQHeader>(daqHandle, daqHandle->size()-1);
  merge::ScanInfo scan;
  MergeScanned(scan, daq); 
  Scan_coll->push_back(scan);
  evt.put(Scan_coll);

 return;
}



void MergeScan::MergeScanned(merge::ScanInfo& scan, edm::Ptr<raw::DAQHeader> daq)
{

  time_t spilltime = daq->GetTimeStamp();//time info. from DAQ480 software
  tm *timeinfo = localtime(&spilltime);
  
  char scanfilename[20];

  sprintf(scanfilename,"/argoneut/app/users/spitz7/larsoft_new/scan_spitz7.00649.txt");
  std::ifstream scanfile(scanfilename);
  
  if(!scanfile.is_open()){
    std::cerr << "MergeScan:  Could not open file named " << scanfilename << std::endl;
    return;
  }

  int isneutrino,isnotneutrino,ismaybeneutrino,trackind,trackcol,vertindtime,vertcoltime,vertindwire,vertcolwire,numshower;
  std::string first,k;
  int event,run;
    
  while(getline(scanfile,k)){
    
    std::istringstream ins3;
    ins3.clear();
    ins3.str(k);
        ins3>>run>>event>>isnotneutrino>>ismaybeneutrino>>isneutrino>>trackind>>trackcol>>vertindtime>>vertcoltime>>vertindwire>>vertcolwire>>numshower;
  
    if((daq->GetRun()==run)&&(daq->GetEvent()==event ) ){ foundscaninfo=true;}    
        
    if(foundscaninfo){      
      scan.SetIsNeutrino(isneutrino);
      scan.SetIsnotNeutrino(isnotneutrino);
      scan.SetIsMaybeNeutrino(ismaybeneutrino);
      scan.SetTrackInd(trackind);
      scan.SetTrackCol(trackcol);      
      scan.SetVertIndTime(vertindtime);
      scan.SetVertColTime(vertcoltime);
      scan.SetVertIndWire(vertindwire);
      scan.SetVertColWire(vertcolwire);
      scan.SetNumShower(numshower);

      break;
    }
  }
  scanfile.close();
  
  return;
  
}




}
