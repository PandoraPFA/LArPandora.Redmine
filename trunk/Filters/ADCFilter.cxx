////////////////////////////////////////////////////////////////////////
//
// ADCFilter class
//
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

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
#include "ADCFilter.h"
#include "RawData/raw.h"
#include "Filters/ChannelFilter.h"

namespace filter {

   //-------------------------------------------------
   ADCFilter::ADCFilter(fhicl::ParameterSet const & pset)  
   {   
      this->reconfigure(pset);
   }

   //-------------------------------------------------
   ADCFilter::~ADCFilter()
   {
   }
  
   //-------------------------------------------------
   void ADCFilter::reconfigure(fhicl::ParameterSet const& p)
   {
      fDigitModuleLabel = p.get< std::string > ("DigitModuleLabel"); 
      fMinADC           = p.get< double      > ("MinADC");         
   } 

   //-------------------------------------------------
   void ADCFilter::beginJob()
   {
   }

   //-------------------------------------------------
   bool ADCFilter::filter(art::Event &evt)
   { 
      //Read in raw data
      art::Handle< std::vector<raw::RawDigit> > rawdigitHandle;
      evt.getByLabel(fDigitModuleLabel,rawdigitHandle);

      //Make a ChannelFilter
      filter::ChannelFilter *chanFilt = new filter::ChannelFilter();

      if(!rawdigitHandle->size()) return false;

      for(unsigned int i = 0; i<rawdigitHandle->size(); ++i){
         art::Ptr<raw::RawDigit> digit(rawdigitHandle,i);
         unsigned int channel = digit->Channel();
         if(chanFilt->BadChannel(channel)) continue;
         
         //get ADC values after decompressing
         std::vector<short> rawadc(digit->Samples());
         raw::Uncompress(digit->fADC,rawadc,digit->Compression());
         short max = *std::max_element(rawadc.begin(),rawadc.end()) - digit->GetPedestal();
         if(max>=fMinADC) return true;//found one ADC value above threshold, pass filter  
      }

      return false;//didn't find ADC above threshold

   }
	      
} //end namespace
