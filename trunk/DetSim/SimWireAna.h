////////////////////////////////////////////////////////////////////////
// $Id: SimWire.h,v 1.12 2010/04/23 20:30:53 seligman Exp $
//
// SimWire class designed to simulate signal on a wire in the TPC
//
// echurch@fnal.gov, brebel@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////
#ifndef SIMWIREANA_H
#define SIMWIREANA_H

#include "RawData/raw.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArFFT.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH2.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace detsim {

  /// Base class for creation of raw signals on wires. 
  class SimWireAna : public art::EDAnalyzer {
    
  public:
        
    explicit SimWireAna(fhicl::ParameterSet const& pset); 
    virtual ~SimWireAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::string            fDetSimModuleLabel;///< name of module that produced the digits
    TH1F*                  fDiffs;            ///< histogram of Raw tdc to tdc differences

    TH1F*                  fCompressErr;      ///< histogram of difference between original 
                                              ///<tdc value and compressesed value
    TH1F*                  fCompressFactor;   ///< compression factor 

    TH2F*                  fRawVsCompress;    ///< histogram of original tdc value vs compressesed value
    TH2F*                  fCompressErr2D;    ///< histogram of original tdc value vs compressesed value
    

  }; // class SimWire

}

#endif // SIMWIREANA_H
