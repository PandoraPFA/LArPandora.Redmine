////////////////////////////////////////////////////////////////////////
/// 
/// \brief Module to characterize KINEMATICS
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "Utilities/LArProperties.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <TH1F.h>
#include <TH2F.h>
class TTree;
namespace kin {
   
 class Kinematics :  public art::EDAnalyzer {
    
  public:
    
    explicit Kinematics(fhicl::ParameterSet const& pset); 
    virtual ~Kinematics();        

    void analyze (const art::Event& evt);
    void beginJob();
    
  private:
    TTree* ftree;
    TH2F*  flepangle;
    std::string fDBScanModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fVertexModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fScanModuleLabel;
    std::string fMINOSModuleLabel;
    std::string fLarTracks_label;
    int fm_event, fm_pdgcode;
    float fm_leptheta_reco,fm_leptheta_true,fm_lepphi_reco,fm_lepphi_true,fm_init_x_reco,fm_init_y_reco,fm_init_z_reco,fm_init_x_true,fm_init_y_true,fm_init_z_true,fm_tpcexit_x_true,fm_tpcexit_y_true,fm_tpcexit_z_true,fm_tpcexit_x_reco,fm_tpcexit_y_reco,fm_tpcexit_z_reco;
  };
    
}



#endif // KINEMATICS_H
