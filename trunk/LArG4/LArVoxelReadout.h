////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.h
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \version $Id: LArVoxelReadout.h,v 1.2 2009/03/31 17:58:39 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// One way to implement voxels in Geant4 is to create a parallel
/// "read-out" geometry along with the real, physical geometry.  The
/// read-out geometry is implemented in LArVoxelReadoutGeometry; this
/// class is the sensitive detector for that geometry.  That is,
/// Geant4 will call this routine every time there is a step within a
/// volume of the read-out geometry; this routine then accumulates
/// information from that step.
///
/// In general, Geant4 expects to have per-event user information
/// attached to the G4Event in some way; their G4VSensitiveDetector
/// class supports this by allowing user-defined hit collections to
/// added to a G4HCOfThisEvent object (a collection of hit
/// collections; yes, it makes my head ache too!) that's part of each
/// G4Event.  
///
/// This class works differently, by accumulating the information in
/// its internal sim::LArVoxelList.  See LArVoxelListAction for how
/// this information is made available to the main LArG4 module.
/// 
/// Why define a parallel geometry?  Here are some reasons:
///
/// - The regular LAr TPC is one large volume of liquid argon.  When
///   Geant4 does its physics modeling, it can be unconstrained in
///   step size by the voxels.  Only for readout would the steps be
///   sub-divided.
///
/// - There may be more than one kind of readout, depending on a
///   detector's instrumentation (e.g., PMTs in addition to the wire
///   planes).  It's possible that the voxelization appropriate for
///   the wire planes may not be an appropriate readout for the other
///   readouts.  Geant4 allows the construction of multiple parallel
///   readouts, so this mechanism is relatively easy to extend for
///   each type of readout.

#ifndef LArG4_LArVoxelReadout_h
#define LArG4_LArVoxelReadout_h

#include <G4VSensitiveDetector.hh>
#include <globals.hh>

#include "Simulation/LArVoxelList.h"
#include "Simulation/SimChannel.h"

// Forward declarations
class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

namespace larg4 {

  class LArVoxelReadout : public G4VSensitiveDetector
  {
  public:
    // Constructor.  Note that the typical argument for a sensitive
    // detector is the name; it's omitted here, but assigned in
    // LArVoxelReadout.cc.
    explicit LArVoxelReadout();

    // Destructor
    virtual ~LArVoxelReadout();

    // Required for classes that inherit from G4VSensitiveDetector.
    //
    // Called at start and end of each event.
    virtual void Initialize(G4HCofThisEvent*);
    virtual void EndOfEvent(G4HCofThisEvent*);

    // Called to clear any accumulated information.
    virtual void clear();

    // The key method of this class.  It's called by Geant4 for each
    // step within the read-out geometry.  It accumulates the energy
    // in the G4Step in the LArVoxelList.
    virtual G4bool ProcessHits( G4Step*, G4TouchableHistory* );

    // Empty methods; they have to be defined, but they're rarely
    // used in Geant4 applications.
    virtual void DrawAll();
    virtual void PrintAll();

    // Independent method; returns the accumulated information
    const std::vector<sim::SimChannel>& GetSimChannels() const { return fChannels; }

  private:

    void DriftIonizationElectrons(double energy, 
				  double dx,
				  G4ThreeVector stepMidPoint,
				  int trackID);

    // Used in electron-cluster calculations
    // External parameters for the electron-cluster calculation.
    // obtained from LArG4Parameters, LArProperties, and DetectorProperties services
    double fDriftVelocity;	       
    double fRecombA;	       
    double fRecombk;	       
    double fLongitudinalDiffusion;
    double fTransverseDiffusion;  
    double fElectronLifetime;     
    double fElectronClusterSize;  
    double fGeVToElectrons;       
    double fSampleRate; 	       
    int    fTriggerOffset;        

    std::map<unsigned int, sim::SimChannel  > fChannelMap; ///< Map of channel number to SimChannel object
    std::vector<sim::SimChannel>            fChannels;     ///< Collection of SimChannels
  };

}

#endif // LArG4_LArVoxelReadout_h
