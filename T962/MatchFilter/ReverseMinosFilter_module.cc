////////////////////////////////////////////////////////////////////////
//
// ReverseMInosFilter_module class:
// Algoritm to produce a filtered event file having
// events with at least one MINOS track point back to T962 TPC.
//
//
////////////////////////////////////////////////////////////////////////
#ifndef REVERSEMINOSFILTER_H
#define REVERSEMINOSFILTER_H

#include <string>

/// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
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
#include "T962/T962_Objects/MINOS.h"



namespace filt {

  class ReverseMinosFilter : public art::EDFilter  {
    
  public:
    
    explicit ReverseMinosFilter(fhicl::ParameterSet const& pset); 
    virtual ~ReverseMinosFilter();
         
    
    bool filter(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
   

  private:
    
     std::string fMinosTracks_label;

    double fdZ; //max. z distance into MINOS for matchable tracks.

  protected: 

    
  }; // class ReverseMinosFilter

  //-------------------------------------------------
  ReverseMinosFilter::ReverseMinosFilter(fhicl::ParameterSet const & pset)  
  {   
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  ReverseMinosFilter::~ReverseMinosFilter()
  {
  }
  
  //-------------------------------------------------
  void ReverseMinosFilter::reconfigure(fhicl::ParameterSet const& pset)
  {
    fMinosTracks_label = pset.get< std::string >("minostracks");
    fdZ = pset.get<double>("dZ");
  } 

  //-------------------------------------------------
  void ReverseMinosFilter::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    //fSelectedEvents = tfs->make<TH1D>("fSelectedEvents", "Number of Selected Events", 3, 0, 3); //counts the number of selected events 
  }

  //-------------------------------------------------
  bool ReverseMinosFilter::filter(art::Event &evt)
  { 

    art::Handle< std::vector<t962::MINOS> > MinosHandle;
    evt.getByLabel(fMinosTracks_label,MinosHandle);
    if(MinosHandle.failedToGet()){
      mf::LogWarning("filter") << "No MINOS information found with label = "
                                 << fMinosTracks_label << ". Skipping.\n";
      return false;
    }
    
    double D=(90*0.5)+(42.4*2.54)-5.588; //147.108 cm is distance from the front (upstream) of the TPC to the 1st Minos plane
    double x_offset = 117.4;
    double y_offset = 19.3;
    
    if(MinosHandle->size()>0){
      for(unsigned int i = 0; i<MinosHandle->size();++i){
        art::Ptr<t962::MINOS> minos_track(MinosHandle,i);
        if((100.0*minos_track->ftrkVtxZ)>fdZ) continue;//MINOS track too far into detector
        double x_start = 100.0*minos_track->ftrkVtxX - x_offset;//in T962 coordinates
        double y_start = 100.0*minos_track->ftrkVtxY + y_offset;//in T962 coordinates
        double z_start = D+(100.0*minos_track->ftrkVtxZ);//in T962 coordinates

        //walk track backwards along its directional cosines towards T962
        for(unsigned int j = 0; j<500; ++j){
          //move backwards in 1cm increments along track direction
          double x = x_start - j*minos_track->ftrkdcosx;
          double y = y_start - j*minos_track->ftrkdcosy;
          double z = z_start - j*minos_track->ftrkdcosz;
          if(x>0.0 && x<47.5 && y>-20.0 && y<20.0 && z>0.0 && z<90.0){
            //std::cout << "MINOS ID#" << minos_track->ftrkIndex << " is a Match Candidate!" << std::endl;
            return true;
          }
          if(z<-5.0) break;
        }
      }
    }
    else{
      return false;
    }

    return false;
     
  } // bool

  DEFINE_ART_MODULE(ReverseMinosFilter);

} //namespace filt

#endif // REVERSEMINOSFILTER_H
