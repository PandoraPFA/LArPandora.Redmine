////////////////////////////////////////////////////////////////////////
//
// MatchFilter class:
// Algoritm to produce a filtered event file having
// events with atleast 1 matched track with MINOS
//
// saima@ksu.edu
//
////////////////////////////////////////////////////////////////////////

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
#include "MatchFilter.h"
#include "T962/T962_Objects/MINOSTrackMatch.h"

namespace filt{

  //-------------------------------------------------
  MatchFilter::MatchFilter(fhicl::ParameterSet const & pset)  
  {   
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  MatchFilter::~MatchFilter()
  {
  }
  
  //-------------------------------------------------
  void MatchFilter::reconfigure(fhicl::ParameterSet const& p)
  {
    fTrackMatchModuleLabel  =  p.get< std::string >("TrackMatchModuleLabel");
  } 

  //-------------------------------------------------
  void MatchFilter::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fSelectedEvents = tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3); //counts the number of selected events 
  }

  //-------------------------------------------------
  bool MatchFilter::filter(art::Event &evt)
  { 
    
  art::Handle< std::vector<t962::MINOSTrackMatch> > trackmatchListHandle;
  evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle);

  art::PtrVector<t962::MINOSTrackMatch> trackmatchlist;
  if(evt.getByLabel(fTrackMatchModuleLabel,trackmatchListHandle))
  for (unsigned int i = 0; i < trackmatchListHandle->size(); i++){
    art::Ptr<t962::MINOSTrackMatch> trackmatchHolder(trackmatchListHandle,i);
    trackmatchlist.push_back(trackmatchHolder);
  }
  
  int nmatched_reco = -99999;
  nmatched_reco = trackmatchlist.size();

  if(nmatched_reco>0){
    fSelectedEvents->Fill(1); 
    //std::cout << "event passed the filter " << evt.id().event() << std::endl;
    return true;
  }
  
  else {
    //std::cout << "event failed the filter " << evt.id().event() << std::endl;
    return false;
  }

  } // bool  
} // namespace
    
//------------------------------------------------   

