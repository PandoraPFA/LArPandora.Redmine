////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to characterize KINEMATICS
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "Utilities/LArProperties.h"

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>

namespace kin {
   
 class Kinematics :  public art::EDProducer {
    
  public:
    
    explicit Kinematics(fhicl::ParameterSet const& pset); 
    virtual ~Kinematics();        

    void produce(art::Event& evt);
    void beginJob();
    
  private:

    TH2F*  findcol;
    std::string fDBScanModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fVertexModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fScanModuleLabel;
    std::string fMINOSModuleLabel;
    
  };
    
}



#endif // KINEMATICS_H
