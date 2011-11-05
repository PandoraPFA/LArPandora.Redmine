
////////////////////////////////////////////////////////////////////////
// $Id: CalWireAna.h,v 1.12 2010/04/23 20:30:53 seligman Exp $
//
// CalWireAna class designed to make histograms
//
// echurch@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////
#ifndef CALWIREANA_H
#define CALWIREANA_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <string>

class TH2F;
class TH1F;

namespace geo { class Geometry; }

namespace caldata {

  /// Base class for creation of raw signals on wires. 
  class CalWireAna : public art::EDAnalyzer {
    
  public:
        
    explicit CalWireAna(fhicl::ParameterSet const& pset); 
    virtual ~CalWireAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void endJob();

  private:

    std::string            fCalWireModuleLabel;///< name of module that produced the wires
    std::string            fDetSimModuleLabel; //< name of module that produced the digits

    TH1F*                  fDiffsR;            ///< histogram of Raw tdc to tdc differences
    TH1F*                  fDiffsW;            ///< histogram of Wire (post-deconvoution) tdc to tdc differences
    TH1F*                  fDiffsRW; 
    TH1F*                  fDiffsRWgph; 
    TH1F*                  fMin; 
    TH1F*                  fMax; 
    TH1F*                  fIR; 
    TH1F*                  fCR; 
    TH1F*                  fIW; 
    TH1F*                  fCW; 
    TH1F*                  fRawIndPeak;      //Raw Induction Peak Values
    TH1F*                  fRawColPeak;      //Raw Collection Peak Values
    TH1F*                  fCalIndPeak;      //Calibrated Induction Peak Values
    TH1F*                  fCalColPeak;      //Calibrated Collection Peak Values
    TH1F*                  fNoiseHist;       //Noise Frequency Spectrum
    TH1F*                  fNoiseRMS;        //Noise RMS values
    TH2F*                  fWireSig;
    TH2F*                  fRawSig;
    TH2F*                  fRD_WireMeanDiff2D;    ///< histogram of difference between original tdc value and compressesed value vs original value
    TH2F*                  fRD_WireRMSDiff2D;    ///< histogram of difference between original tdc value and compressesed value vs original value
    TH2F*                  fDiffsRWvsR;
    TH2F*                  fDiffsRWvsRgph;
    TH2F*                  fWindow;
    

  }; // class CalWireAna

} // End caldata namespace.

#endif // CALWIREANA_H
