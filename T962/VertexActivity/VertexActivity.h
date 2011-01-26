////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to characterize vertex activity
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef VertexActivity_H
#define VertexActivity_H

#include "Utilities/LArProperties.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>

namespace vertex {
   
 class VertexActivity :  public edm::EDProducer {
    
  public:
    
    explicit VertexActivity(edm::ParameterSet const& pset); 
    ~VertexActivity();        

    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
    
  private:

    TH2F*  fEfficiency;
    TH2F*  fIndEfficiency;
    TH2F*  fColEfficiency;
    TH2F*  findcol;
    std::string fDBScanModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fVertexModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fScanModuleLabel;
    double fCathodetimelocation, fDelta_Cathodetimelocation;
    double fE_lifetime, fDelta_E_lifetime;
    double fRecombination_factor, fDelta_Recombination_factor;
    double fWorkfunction_factor, fDelta_Workfunction_factor;
    double fCalibration_factor, fDelta_Calibration_factor;
    double fActivityRadius;//units are 'ticks'
    
  };
    
}



#endif // VertexActivity_H
