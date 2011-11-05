////////////////////////////////////////////////////////////////////////
/// \file  LArG4Ana.h
/// \brief Check of Geant4 to run the LArSoft detector simulation
///
/// \version $Id: LArG4.h,v 1.11 2010/06/04 21:47:27 bjpjones Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARG4_MNSP_H
#define LARG4_MNSP_H 

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4WrapperProcess.hh"

#include "G4MuNuclearInteraction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4DataQuestionaire.hh"

namespace larg4 {

class MuNuclearSplittingProcess : public G4WrapperProcess {
// Override PostStepDoIt method 
  public:
    MuNuclearSplittingProcess() {};
    ~MuNuclearSplittingProcess() {};
    
    void SetNSplit(G4int nTrx) {fNSplit = nTrx;};
    void SetIsActive(G4bool doIt) {fActive = doIt;}; 
    
  private:
// Data members
    G4int fNSplit;
    G4bool fActive; 
    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);

}; 


}// end namespace

#endif // MNSP
