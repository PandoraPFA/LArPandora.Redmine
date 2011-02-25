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

///Geant4 interface 
namespace T962G4 {  
 
  class T962G4Ana : public art::EDAnalyzer{
  public:

    explicit T962G4Ana(fhicl::ParameterSet const& pset);
    virtual ~T962G4Ana();
    void analyze (const art::Event& evt); 
    void beginJob();

  private:

    std::string fG4ModuleLabel;  // name of the process module label that produced the input 
    
  };

} // namespace T962G4

#endif // T962G4_T962G4_H
