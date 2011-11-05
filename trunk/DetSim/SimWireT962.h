////////////////////////////////////////////////////////////////////////
// $Id: SimWireT962.h,v 1.12 2010/04/23 20:30:53 seligman Exp $
//
// SimWireT962 class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////
#ifndef SIMWIRET962_H
#define SIMWIRET962_H

#include "RawData/raw.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArFFT.h"
#include "Utilities/LArProperties.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"

#include "art/Framework/Core/EDProducer.h"

#include <vector>
#include <string>

namespace art {
  class Event;
  class ParameterSet;
}

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires. 
  class SimWireT962 : public art::EDProducer {
    
  public:
        
    explicit SimWireT962(fhicl::ParameterSet const& pset); 
    virtual ~SimWireT962();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void         ConvoluteResponseFunctions(); ///< convolute electronics and field response
    
    void         SetFieldResponse();           ///< response of wires to field
    void         SetElectResponse();           ///< response of electronics
    
    void         GenNoise(std::vector<double>& array);

    bool                   fResponseSet;      ///< flag of whether to set the response functions or not
    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    std::string            fResponseFile;     ///< response file for induction planes
    raw::Compress_t        fCompression;      ///< compression type to use

    double                 fNoiseFact;        ///< noise scale factor 
    double                 fNoiseWidth;       ///< exponential noise width (kHz) 
    double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
    int                    fNTicks;           ///< number of ticks of the clock
    int                    fNFieldBins;       ///< number of bins for field response
    double                 fSampleRate;       ///< sampling rate in ns
    double                 fCol3DCorrection;  ///< correction factor to account for 3D path of 
                                              ///< electrons thru wires
    double                 fInd3DCorrection;  ///< correction factor to account for 3D path of 
                                              ///< electrons thru wires
    double                 fColFieldRespAmp;  ///< amplitude of response to field 
    double                 fIndFieldRespAmp;  ///< amplitude of response to field 
    std::vector<double>    fShapeTimeConst;   ///< time constants for exponential shaping
    int                    fTriggerOffset;    ///< (units of ticks) time of expected neutrino event
    unsigned int           fNElectResp;       ///< number of entries from response to use
    
    std::vector<double>    fColFieldResponse; ///< response function for the field @ collection plane
    std::vector<double>    fIndFieldResponse; ///< response function for the field @ induction plane
    std::vector<TComplex>  fColShape;         ///< response function for the field @ collection plane
    std::vector<TComplex>  fIndShape;         ///< response function for the field @ induction plane
    std::vector<double>    fChargeWork;
    std::vector<double>    fElectResponse;    ///< response function for the electronics
    std::vector< std::vector<double> > fNoise;///< noise on each channel for each time
    
    TH1D*                fIndFieldResp;     ///< response function for the field @ induction plane
    TH1D*                fColFieldResp;     ///< response function for the field @ collection plane
    TH1D*                fElectResp;        ///< response function for the electronics
    TH1D*                fColTimeShape;     ///< convoluted shape for field x electronics @ col plane
    TH1D*                fIndTimeShape;     ///< convoluted shape for field x electronics @ ind plane
    TH1D*                fNoiseDist;        ///< distribution of noise counts

  }; // class SimWireT962

}

#endif // SIMWIRET962_H
