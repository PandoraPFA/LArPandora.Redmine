////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadout.cxx
/// \brief A Geant4 sensitive detector that accumulates voxel information.
///
/// \version $Id: LArVoxelReadout.cxx,v 1.4 2009/08/12 20:52:14 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "LArG4/LArVoxelReadout.h"
#include "LArG4/ParticleListAction.h"
#include "Geometry/geo.h"
#include "Simulation/LArG4Parameters.h"
#include "Simulation/Particle.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <G4HCofThisEvent.hh>
#include <G4TouchableHistory.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4ThreeVector.hh>
#include <Randomize.hh>

#include <iostream>
#include <ctime>
#include <vector> 

namespace larg4 {

  // Constructor.  Note that we force the name of this sensitive
  // detector to be the value expected by LArVoxelListAction.
  LArVoxelReadout::LArVoxelReadout()
    : G4VSensitiveDetector("LArVoxelSD")
  {
    // Initialize values for the electron-cluster calculation. 
    fChannelMap.clear();
    fChannels.clear();

    // Fetch all parameters in advance. (This eliminates the possilibity of 
    // changing them interactively via the event display, but it's
    // unlikely that the Monte Carlo would be re-run from there.)
    art::ServiceHandle<util::LArProperties>  larp;

    fElectronLifetime = larp->ElectronLifetime();
    fDriftVelocity    = larp->DriftVelocity(larp->Efield(),larp->Temperature())/1000.;

    art::ServiceHandle<sim::LArG4Parameters> lgp;
    fElectronClusterSize   = lgp->ElectronClusterSize();
    fLongitudinalDiffusion = lgp->LongitudinalDiffusion();
    fTransverseDiffusion   = lgp->TransverseDiffusion();
    fGeVToElectrons        = lgp->GeVToElectrons();
    fRecombA               = lgp->RecombA();
    fRecombk 		   = lgp->Recombk();

    art::ServiceHandle<util::DetectorProperties> detprop;
    fSampleRate       = detprop->SamplingRate();
    fTriggerOffset    = detprop->TriggerOffset();

  }

  //-------------------------------------------------------
  LArVoxelReadout::~LArVoxelReadout() {}

  // Called at the start of each event.
  void LArVoxelReadout::Initialize(G4HCofThisEvent*)
  {
    fChannelMap.clear();
    fChannels.clear();
  }

  // Called at the end of each event.
  void LArVoxelReadout::EndOfEvent(G4HCofThisEvent*)
  {
    // iterate over the map and fill the vector of channels
    std::map<unsigned int, sim::SimChannel>::const_iterator itr = fChannelMap.begin();
    while( itr != fChannelMap.end() ){

      fChannels.push_back((itr->second));
      itr++;
    }

  }

  //-----------------------------------------------------------------
  void LArVoxelReadout::clear()
  {
    fChannelMap.clear();
    fChannels.clear();
  }

  // Called for each step.
  G4bool LArVoxelReadout::ProcessHits( G4Step* step, G4TouchableHistory* )
  {
    // All work done for the "parallel world" "box of voxels" in
    // LArVoxelReadoutGeometry makes this a fairly simple routine.
    // First, the usual check for non-zero energy:

    G4double energy = step->GetTotalEnergyDeposit();
    if ( energy > 0 ){
      // The step can be no bigger than the size of the voxel,
      // because of the geometry set up in LArVoxelGeometry and the
      // transportation set up in PhysicsList.  Find the mid-point
      // of the step.
      G4ThreeVector midPoint = 0.5*( step->GetPreStepPoint()->GetPosition() 
				     + step->GetPostStepPoint()->GetPosition() );
      
      // Find the Geant4 track ID for the particle responsible for depositing the
      // energy.  if we are only storing primary EM shower particles, and this energy
      // is from a secondary etc EM shower particle, the ID returned is the primary
      const int trackID = ParticleListAction::GetCurrentTrackID();
      
      // Now calculate the number of ionization electrons produced in this step
      // and drift them to the readout wires
      G4ThreeVector totstep = step->GetPostStepPoint()->GetPosition() - step->GetPreStepPoint()->GetPosition();

      // Note that if there is no particle ID for this energy deposit, the
      // trackID will be sim::NoParticleId.
      DriftIonizationElectrons(energy/MeV, totstep.mag()/cm, midPoint, trackID);

    } // end if energy > 0

    return true;
  }

  //----------------------------------------------------------------------------
  // energy is passed in with units of MeV, dx has units of cm
  void LArVoxelReadout::DriftIonizationElectrons(double energy, 
						 double dx,
						 G4ThreeVector stepMidPoint,
						 int trackID)
  {

    // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
    // R = A/(1 + (dE/dx)*k)
    // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
    // from GeV/voxel width
    // A = 0.800 +/- 0.003
    // k = (0.097+/-0.001) g/(MeVcm^2)
    
    // the dx depends on the trajectory of the step
    
    // This routine gets called frequently, once per every particle
    // traveling through every voxel. Use whatever tricks we can to
    // increase its execution speed.

    static double LifetimeCorr_const = -1000. * fElectronLifetime;
    static double nElectrons_const   = fGeVToElectrons  *1.e-3;
    static double LDiff_const        = sqrt(2.*fLongitudinalDiffusion);
    static double TDiff_const        = sqrt(2.*fTransverseDiffusion);
    static double RecipDriftVel      = 1./fDriftVelocity;
    
    // Map of electrons to store - catalogued by map[channel][tdc]
    static std::map<unsigned int, std::map<unsigned int,double> >  ElectronsToStore;

    double recomb = fRecombA/(1. + (energy/dx)*fRecombk);
        
    static double xyz1[3] = {0.};

    double xyz[3] = {stepMidPoint.x() / cm,
		     stepMidPoint.y() / cm,
		     stepMidPoint.z() / cm};

    // figure out which TPC we are in based on the stepMidPoint
    art::ServiceHandle<geo::Geometry> geom;
    try{
      const geo::TPCGeo &tpc = geom->PositionToTPC(xyz);

      // X drift distance - the drift direction can be either in 
      // the positive or negative direction, so use fabs
      // \todo fix this up for use with Bo as it drifts in the y direction
      double XDrift = fabs(stepMidPoint.x()/cm - tpc.PlaneLocation(0)[0]);
      
      if(XDrift < 0.) return; 
      
      // Drift time (nano-sec)
      double TDrift = XDrift * RecipDriftVel;
      double lifetimecorrection = TMath::Exp(TDrift/LifetimeCorr_const);
      double nElectrons = lifetimecorrection * energy * recomb * nElectrons_const;
      // Longitudinal & transverse diffusion sigma (cm)
      double SqrtT = sqrt(TDrift);
      double LDiffSig = SqrtT * LDiff_const;
      double TDiffSig = SqrtT * TDiff_const;
      int nClus = 1 + (int)(nElectrons/fElectronClusterSize);
      
      // Compute arrays of values as quickly as possible.
      std::vector< double > XDiff(nClus); 
      std::vector< double > YDiff(nClus);
      std::vector< double > ZDiff(nClus);
      std::vector< double > nElDiff(nClus,fElectronClusterSize);
      // Drift the last cluster with smaller size
      nElDiff[nClus-1] = nElectrons - (nClus-1)*fElectronClusterSize;
      // Smear drift times by x position and drift time
      G4RandGauss::shootArray( nClus, &XDiff[0], 0., LDiffSig);
      // Smear the Y,Z position by the transverse diffusion
      G4RandGauss::shootArray( nClus, &YDiff[0], stepMidPoint.y()/cm,TDiffSig);
      G4RandGauss::shootArray( nClus, &ZDiff[0], stepMidPoint.z()/cm,TDiffSig);
      
      // make a collection of electrons for each plane
      for(size_t p = 0; p < tpc.Nplanes(); ++p){
	
	double Plane0Pitch = tpc.Plane0Pitch(p);
	
	xyz1[0] = tpc.PlaneLocation(0)[0] - Plane0Pitch;
	//changed to "-" sign due to the change of Plane0Pitch output to positive. Andrzej
	
	// Drift nClus electron clusters to the induction plane
	for(int k = 0; k<nClus; ++k){
	  // Correct drift time for longitudinal diffusion and plane
	  double TDiff = TDrift + (Plane0Pitch + XDiff[k]) * RecipDriftVel;
	  xyz1[1] = YDiff[k];
	  xyz1[2] = ZDiff[k];

	  // \todo think about effects of drift between planes
	  
	  // grab the nearest channel to the xyz position
	  try{
	    unsigned int channel = geom->NearestChannelFast(xyz1, p, 0);
	    
	    // \todo check on what happens if we allow the tdc value to be 
	    // beyond the end of the expected number of ticks

	    unsigned int tdc = (unsigned int)(TDiff/fSampleRate) + fTriggerOffset;         
	   
	    // Add electrons produced by each cluster to the map
	    ElectronsToStore[channel][tdc]+=nElDiff[k];
	  }
	  catch(cet::exception &e){
	    mf::LogWarning("LArVoxelReadout") << "unable to drift electrons from point (" 
					      << xyz[0] << "," << xyz[1] << "," << xyz[2]
					      << ") with exception " << e; 
	  }	  
	} // end loop over clusters	
      } // end loop over planes

      // Now store them in SimChannels

      // check if the current channel is already in the map, otherwise add it
      for(std::map<unsigned int,std::map<unsigned int,double> >::const_iterator itread = ElectronsToStore.begin();
	  itread!=ElectronsToStore.end();
	  itread++){
	unsigned int channel = itread->first;
	std::map<unsigned int, sim::SimChannel>::iterator itchannelmap = fChannelMap.find(channel);
	if( itchannelmap != fChannelMap.end() ) 
	  for(std::map<unsigned int, double>::const_iterator itreadinner = itread->second.begin();
	      itreadinner!=itread->second.end();
	      itreadinner++)
	    itchannelmap->second.AddIonizationElectrons(trackID, itreadinner->first, itreadinner->second, xyz);
	else{
	  sim::SimChannel sc(channel);
	  for(std::map<unsigned int, double>::const_iterator itreadinner = itread->second.begin();
	      itreadinner!=itread->second.end();
	      itreadinner++)
	    sc.AddIonizationElectrons(trackID, itreadinner->first, itreadinner->second, xyz);
	  fChannelMap[channel] = sc;
	}
      }
      ElectronsToStore.clear();
      
    } // end try intended to catch points where TPC can't be found
    catch(cet::exception &e){
      mf::LogWarning("LArVoxelReadout") << "step cannot be found in a TPC\n"
					<< e;
    }

    return;
  }

  // Never used.
  void LArVoxelReadout::DrawAll()  {}
  void LArVoxelReadout::PrintAll() {}

} // namespace larg4
