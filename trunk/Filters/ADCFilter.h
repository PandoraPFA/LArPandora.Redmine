////////////////////////////////////////////////////////////////////////
//
// ADCFilter class:
// Algorithm to ignore events with no ADC values 
// above user-defined threshold.
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef ADCFILTER_H
#define ADCFILTER_H

#include "art/Framework/Core/EDFilter.h"

namespace filter {

   class ADCFilter : public art::EDFilter  {
    
   public:
    
      explicit ADCFilter(fhicl::ParameterSet const& ); 
      virtual ~ADCFilter();
         
    
      bool filter(art::Event& evt);
      void reconfigure(fhicl::ParameterSet const& p);
      void beginJob();
   

   private: 
 
      std::string fDigitModuleLabel;
      double      fMinADC;  
 
   protected: 
    
   }; // class ADCFilter

}

#endif // ADCFILTER_H
