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
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Core/FindMany.h"


//Larsoft Includes
#include "T962/MatchFilter/MatchFilter.h"
#include "T962/T962_Objects/classes.h"
#include "RecoBase/recobase.h"

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
     fTracks_label           = p.get< std::string >("LArTracksModuleLabel");
     fTrackMatchModuleLabel  = p.get< std::string >("TrackMatchModuleLabel");
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

     art::Handle< std::vector<recob::Track> > LarTrackHandle;
     evt.getByLabel(fTracks_label,LarTrackHandle);

     art::FindOne<t962::MINOS> fomatch(LarTrackHandle, evt, fTrackMatchModuleLabel);

     if(fomatch.size()>0){
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

