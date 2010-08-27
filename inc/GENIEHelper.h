////////////////////////////////////////////////////////////////////////
/// \file  GENIEHelper.h
/// \brief Wrapper for generating neutrino interactions with GENIE
///
/// \version $Id: GENIEHelper.h,v 1.6 2010/04/27 19:03:42 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVGB_GENIEHELPER_H
#define EVGB_GENIEHELPER_H

#include <vector>
#include <set>
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GeomAnalyzerI.h"
#include "EVGDrivers/GMCJDriver.h"

class TH1D;

namespace edm {
  class ParameterSet;
}

namespace simb  { class MCTruth;     }
namespace simb  { class MCFlux;      }
namespace genie { class EventRecord; }

namespace evgb{

  class GENIEHelper {
    
  public:
  
    explicit GENIEHelper(edm::ParameterSet const& pset);
    ~GENIEHelper();

    void                   Initialize();
    bool                   Stop();
    bool                   Sample(simb::MCTruth &truth, 
				  simb::MCFlux  &flux);
    double                 TotalHistFlux();
    double                 POTUsed()          const { return fPOTUsed;     }
    std::string            FluxType()         const { return fFluxType;    }
    std::string            DetectorLocation() const { return fDetLocation; }

  private:

    void InitializeGeometry();
    void InitializeFluxDriver();
    void PackNuMIFlux(simb::MCFlux &flux);
    void PackSimpleFlux(simb::MCFlux &flux);
    void PackMCTruth(genie::EventRecord *record, simb::MCTruth &truth);

    genie::GeomAnalyzerI*    fGeomD;       
    genie::GFluxI*           fFluxD;
    genie::GMCJDriver*       fDriver;

    std::string              fFluxType;          ///histogram or ntuple
    std::string              fFluxFile;          ///name of file containing histograms or ntuples, or txt
    std::string              fBeamName;          ///name of the beam we are simulating
    std::string              fTopVolume;         ///top volume in the ROOT geometry in which to generate events
    std::string              fDetLocation;       ///name of flux window location
    std::vector<TH1D *>      fFluxHistograms;    ///histograms for each nu species

    double                   fTargetA;           ///A of the target nucleus
    double                   fEventsPerSpill;    ///number of events to generate in each spill if not using POT/spill
    double                   fPOTPerSpill;       ///number of pot per spill
    double                   fHistEventsPerSpill;///number of events per spill for histogram fluxes - changes each spill
    double                   fSpillTotal;        ///total of either pot or events for this spill
    double                   fMonoEnergy;        ///energy of monoenergetic neutrinos
    double                   fPOTUsed;           ///pot used from flux ntuple
    double                   fXSecMassPOT;       ///product of cross section, mass and POT/spill for histogram fluxes
    double                   fTotalHistFlux;     ///total flux of neutrinos from flux histograms for used flavors
    TVector3                 fBeamDirection;     ///direction of the beam for histogram fluxes
    TVector3                 fBeamCenter;        ///center of beam for histogram fluxes - must be in meters
    double                   fBeamRadius;        ///radius of cylindar for histogram fluxes - must be in meters
    double                   fDetLength;         ///length of the detector in meters
    double                   fDetectorMass;      ///mass of the detector in kg
    double                   fSurroundingMass;   ///mass of material surrounding the detector that is intercepted by 
                                                 ///the cylinder for the histogram flux in kg
    double                   fGlobalTimeOffset;  ///overall time shift (ns) added to every particle time
    double                   fRandomTimeOffset;  ///additional random time shift (ns) added to every particle time 
    double                   fZCutOff;           ///distance in z beyond the end of the detector that you allow interactions, in m
    std::set<int>            fGenFlavors;        ///pdg codes for flavors to generate
    std::vector<std::string> fEnvironment;       ///environmental variables and settings used by genie
  };
}
#endif //EVGB_GENIEHELPER_H
