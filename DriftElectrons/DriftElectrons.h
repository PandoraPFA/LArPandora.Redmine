////////////////////////////////////////////////////////////////////////
/// \file  DriftElectrons.h
/// \brief Module to drift ionization electrons to wire planes
///
/// \version $Id: DriftElectrons.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  baller@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef DRIFTELECTRONS_H
#define DRIFTELECTRONS_H

#include "Utilities/LArProperties.h"

#include "art/Framework/Core/EDProducer.h"
#include "TRandom3.h"
#include <vector>
#include <string>

class TH1D;
  
namespace sim {
  class Electrons;
}

///Drifting of electrons
namespace dfe {
  
  ///class describing the drift of electrons in LAr
  class DriftElectrons : public art::EDProducer {
    
  public:
    
    explicit DriftElectrons(fhicl::ParameterSet const& pset);
    virtual ~DriftElectrons();
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet p);

  private:
        
    TRandom3 fRandom;              ///< random number generator

    double fRecombA;               ///< factor for calculating recombination
    double fRecombk;               ///< factor for calculating recombination
    double fGeV2Elect;             ///< Conversion from GeV energy loss to electrons
    double fLongDiff;              ///< Longitudinal diffusion coefficient (cm^2/ns)
    double fTranDiff;              ///< Transverse diffusion coefficient (cm^2/ns)
    double fClusterSize;           ///< Max number of electrons to diffuse as a single cluster
    std::string fLArG4ModuleLabel; ///< module that made the voxels

    TH1D*  fChannels;              ///< channels that have electrons drifted to them
    TH1D*  fDiffuseX;              ///< diffusion in the X direction;
    TH1D*  fDiffuseY;              ///< diffusion in the Y direction;
    TH1D*  fDiffuseZ;              ///< diffusion in the Z direction;
    TH1D*  fNumVoxels;             ///< number of voxels from each event
    TH1D*  fVoxelEnergy;           ///< energy from each voxel
    
  }; // class DriftElectrons
}

#endif // DRIFTELECTRONS_H
