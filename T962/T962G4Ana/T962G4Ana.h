////////////////////////////////////////////////////////////////////////
/// \file  T962G4Ana.h
/// \brief Check of Geant4 to run the LArSoft detector simulation
///
/// \version $Id: T962G4.h,v 1.11 2010/06/04 21:47:27 bjpjones Exp $
/// \author  joshua.spitz@yale.edu
/// \echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef T962G4_T962G4ANA_H
#define T962G4_T962G4ANA_H 

#include "art/Framework/Core/EDAnalyzer.h"

#include <cstring>
#include <TTree.h>
class TH1D;
class TTree;
///Geant4 interface 
namespace T962G4 {  
 
  class T962G4Ana : public art::EDAnalyzer{
  public:
   
    explicit T962G4Ana(fhicl::ParameterSet const& pset);
    virtual ~T962G4Ana();
    void analyze (const art::Event& evt); 
    void beginJob();

  private:
    TTree* ftree;
    std::string fG4ModuleLabel;  // name of the process module label that produced the input 
    int fm_run, fm_event, fm_trackID, fm_pdgcode; 
    float fm_mass, fm_energylost, fm_deflectionangle, fm_tpcexit_px, fm_tpcexit_py, fm_tpcexit_pz, fm_tpcexit_x, fm_tpcexit_y, fm_tpcexit_z, fm_tpcexit_energy, fm_minosenter_px, fm_minosenter_py, fm_minosenter_pz, fm_minosenter_x, fm_minosenter_y, fm_minosenter_z, fm_minosenter_energy; 
  };

} // namespace T962G4

#endif // T962G4_T962G4_H
