///////////////////////////////////////////////////////
//
// EventFilter Class
//
//
//  echurch@fnal.gov
//
///////////////////////////////////////////////////////
#include "Filters/EventFilter.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TMath.h"

//Framework Includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


//Larsoft Includes

void filter::EventFilter::beginJob()
{
}

//-------------------------------------------------
void filter::EventFilter::endJob()
{
}

///////////////////////////////////////////////////////

void filter::EventFilter::reconfigure(fhicl::ParameterSet const& p)
{
  fBadEvents  = p.get < std::vector <unsigned int> >("BadEvents");	       
  fBadRuns    = p.get < std::vector <unsigned int> >("BadRuns");	       
}


filter::EventFilter::EventFilter(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

filter::EventFilter::~EventFilter()
{
}

bool filter::EventFilter::filter(art::Event &evt)
{   
  unsigned int evtNo = (unsigned int) evt.id().event();
  unsigned int runNo = (unsigned int) evt.run();
  std::vector <unsigned int> sobe = SetOfBadEvents();
  std::vector <unsigned int> sobr = SetOfBadRuns();
  if (sobe.size() != sobr.size()) 
    {
      throw cet::exception("EventFilter.cxx: ") << " BadEvent and BadRun list must be same length. Line " <<__LINE__ << ", " << __FILE__ << "\n";
    }

  for (unsigned int ii=0; ii<sobe.size(); ++ii)
    {
      if(sobe.at(ii)==evtNo && sobr.at(ii)==runNo) 
	{
	  mf::LogInfo("EventFilter: ") << "\t\n Skipping run/event " << runNo <<"/"<< evtNo << " by request.\n";
	  return false;
	}
    }
  return true;  
}


