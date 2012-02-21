////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to analyze Muons
///
/// \author  msoderbe@syr.edu
////////////////////////////////////////////////////////////////////////

#ifndef RADARGONANA_H
#define RADARGONANA_H

#include "art/Persistency/Common/Ptr.h"
#include "RecoBase/Hit.h"

#include "Utilities/LArProperties.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>

///T962 radioactive argon analysis code
namespace radargon { 

  class RadArgonAna :  public art::EDAnalyzer {

  public:

    explicit RadArgonAna(fhicl::ParameterSet const& pset); 
    virtual ~RadArgonAna();        

    void analyze (const art::Event& evt);
    void beginJob();

  private:

    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fHitSensitivity;
    std::string fHits_label;
    std::string fClusterModuleLabel;
    TH2S*  fPlaneHits;
    TH1S*  fIndPlaneHits;
    TH1S*  fColPlaneHits;
    //TH1S*  fChannel;
    TH1F*  fCharge;
    TH1F*  fSingleHitColCharge;
    TH1F*  fSingleHitIndCharge;
    TH1F*  fPeakTime;
    TH1F*  fTimeDifference;
    int         ftmatch;

  };  
}



#endif // KINEMATICS_H
