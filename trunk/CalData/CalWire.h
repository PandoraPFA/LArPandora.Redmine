////////////////////////////////////////////////////////////////////////
// $Id: CalWire.h,v 1.16 2010/08/06 23:06:26 bpage Exp $
//
// CalWire class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//
////////////////////////////////////////////////////////////////////////
#ifndef CALWIRE_H
#define CALWIRE_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework

#include "TComplex.h"

///creation of calibrated signals on wires
namespace caldata {

  class CalWire : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWire(fhicl::ParameterSet const& pset); 
    virtual ~CalWire();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
 
  private:
    
    std::string  fResponseFile;      ///< response file containing transformed 
                                     ///< shape histograms and decay constants
    int          fDataSize;          ///< size of raw data on one wire
    int          fExpEndBins;        ///< number of end bins to consider for tail fit    
    int          fPostsample;        ///< number of postsample bins
    std::string  fDigitModuleLabel;  ///< module that made digits

    std::vector<std::vector<TComplex> > fKernelR;      ///< holds transformed induction 
                                                       ///< response function
    std::vector<std::vector<TComplex> > fKernelS;      ///< holds transformed induction 
                                                       ///< response function
    std::vector<double>                 fDecayConstsR; ///< vector holding RC decay 
                                                       ///< constants
    std::vector<double>                 fDecayConstsS; ///< vector holding RC decay 
                                                       ///< constants
    std::vector<int>                    fKernMapR;     ///< map telling which channels  
                                                       ///< have which response functions
    std::vector<int>                    fKernMapS;     ///< map telling which channels 
                                                       ///< have which response functions
  protected: 
    
  }; // class CalWire
}
#endif // CALWIRE_H
