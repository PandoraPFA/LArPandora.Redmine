////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to characterize KINEMATICS
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef SIMKINEMATICS_H
#define SIMKINEMATICS_H

#include "Utilities/LArProperties.h"

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>
class TTree;

///T962 analysis for simulation kinematics
namespace simkin {
   
 class SimKinematics :  public art::EDProducer {
    
  public:
    
    explicit SimKinematics(fhicl::ParameterSet const& pset); 
    virtual ~SimKinematics();        

    void produce(art::Event& evt);
    void beginJob();
    
  private:
    TTree* ftree;

    
    int fm_run;          // Run number
    int fm_event;        // Event number
    int fm_CCNC;
    int fm_PDG;
    int fm_mode;
    int fm_hitnuc;
    
 
    float fm_leppx;
    float fm_leppy;
    float fm_leppz;
    float fm_lepE;
    float fm_nuE;
    float fm_W;
    float fm_qsqr;
    float fm_vertx;
    float fm_verty;
    float fm_vertz;
    float fm_lepphi;// range = 0-2pi
    float fm_leptheta;
    
    std::string fDBScanModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fVertexModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fScanModuleLabel;
    std::string fMINOSModuleLabel;
    
  };
    
}



#endif // SIMKINEMATICS_H
