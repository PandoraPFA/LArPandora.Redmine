////////////////////////////////////////////////////////////////////////
/// \file  DriftElectrons.cxx
/// \brief Module to drift ionization electrons to wire planes
///
/// \version $Id: DriftElectrons.cxx,v 1.20 2010/04/30 16:10:37 brebel Exp $
/// \author  baller@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <iostream>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "DriftElectrons/DriftElectrons.h"
#include "Geometry/geo.h"
#include "Simulation/sim.h"
#include "Simulation/LArVoxelCalculator.h"

#include "TH1.h"
#include "TMath.h"
#include <ctime>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace dfe{

  static bool chan_sort(const sim::Electrons e1, const sim::Electrons e2)
  {
    return e1.Channel() < e2.Channel();
  }

  //-----------------------------------------------------
  DriftElectrons::DriftElectrons(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  maximum allowed seed for RandomNumberGenerator
    // is 900000000, so time(0) returns the number of seconds since 1970 which is 
    // too large, so subtract the upper limit from the returned value. As of
    // June 28, 2011, the resulting number is 0.45 the limit, so we should be
    // good for about another 41 years before having to revisit this statement.
    // \todo Is this solution for setting a seed acceptable when submitting jobs to the grid?
    unsigned int seed = pset.get< unsigned int >("Seed", time(0) - 902000000);

    createEngine(seed);

    produces< std::vector<sim::Electrons> >();

    // Conversion from GeV energy loss to electrons
    // 0.7 = recombination factor for 500 V/cm
    // 23.6 eV per ion pair * 1E9 ev/GeV
    fGeV2Elect = 1.E9/23.6 ;

    // There's a factor of two in the diffusion equation
    fLongDiff *= 2.;
    fTranDiff *= 2.;
  }

  //-----------------------------------------------------
  DriftElectrons::~DriftElectrons()
  {
  }

  //-----------------------------------------------------
  void DriftElectrons::reconfigure(fhicl::ParameterSet p)
  {
    fRecombA          = p.get< double      >("RecombA");    		
    fRecombk          = p.get< double 	   >("Recombk");    	
    fLongDiff         = p.get< double 	   >("LongDiff");   	
    fTranDiff         = p.get< double 	   >("TranDiff");   	
    fClusterSize      = p.get< double 	   >("ClusterSize");	
    fLArG4ModuleLabel = p.get< std::string >("LArG4ModuleLabel");

    return;
  }

  //-----------------------------------------------------
  void DriftElectrons::beginJob()
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    art::ServiceHandle<geo::Geometry> geo;

    fChannels    = tfs->make<TH1D>("channels",  ";channel;# electrons",  geo->Nchannels(), 0, geo->Nchannels());
    fDiffuseX    = tfs->make<TH1D>("diffusex",  ";x diffusion (cm);",    200, -5., 5.  );
    fDiffuseY    = tfs->make<TH1D>("diffusey",  ";y diffusion (cm);",    200, -5., 5.  );
    fDiffuseZ    = tfs->make<TH1D>("diffusez",  ";z diffusion (cm);",    200, -5., 5.  );
    fNumVoxels   = tfs->make<TH1D>("numvoxels", ";voxels in event;",     500,  0., 5000.);
    fVoxelEnergy = tfs->make<TH1D>("voxele",    ";energy (GeV); voxels", 500,  0., 0.01  );
  }

  //-----------------------------------------------------
  void DriftElectrons::produce(art::Event& evt)
  {
    
    //get a collection of electrons
    std::auto_ptr<std::vector<sim::Electrons> > ecol(new std::vector<sim::Electrons>);
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geom;
    
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandGaussQ gauss(engine);

    double xyz[3] = {0.};
    double xyz1[3] = {0.};
    int planes = geom->Nplanes();
  
    ///plane 0 is the first induction plane, every plane has a wire 0
    geom->Plane(0).Wire(0).GetCenter(xyz);

    std::vector<double> pitch(planes, 0.);
    for(int p = 1; p < planes; ++p){
      geom->Plane(p).Wire(0).GetCenter(xyz1);
      pitch[p] = pitch[p-1] + xyz1[0]-xyz[0];
      xyz[0] = xyz1[0];
//       std::cout << "pitch " << p << " is " << pitch[p] << std::endl;
    }
    // reset the xyz array to be at the first induction plane
    geom->Plane(0).Wire(0).GetCenter(xyz);

//     std::cout << "In Drift: X =" << xyz[0] << std::endl;
  
    // Read in the LArVoxelList object(s).
    art::Handle< std::vector<sim::LArVoxelData> > vxlistHandle;
    evt.getByLabel(fLArG4ModuleLabel,vxlistHandle);

    fNumVoxels->Fill(vxlistHandle->size());


    double electronlifetime=larp->ElectronLifetime();
    double driftvelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature())/1000.;
    // There's probably only one LArVoxelList per event, but FMWK
    // always reads a vector of pointers.  For each LArVoxelList:
    for(unsigned int i = 0; i < vxlistHandle->size(); ++i){

      // Get the reference to the LArVoxelID in the LArVoxelList.
      art::Ptr<sim::LArVoxelData> voxel(vxlistHandle, i);
      
      // Get the energy calculated by the LArVoxelList method.
      double Energy = voxel->Energy();
      fVoxelEnergy->Fill(Energy);

      // Get the recombination factor for this voxel - Nucl.Instrum.Meth.A523:275-286,2004
      // R = A/(1 + (dE/dx)*k)
      // dE/dx is given by the voxel energy deposition, but have to convert it to MeV/cm
      // from GeV/voxel width
      // A = 0.800 +/- 0.003
      // k = (0.097+/-0.001) g/(MeVcm^2)

      art::ServiceHandle<sim::LArVoxelCalculator> lvc;
      double recomb = fRecombA/(1. + Energy*(1e3/lvc->VoxelSizeX())*fRecombk);

      
      // X drift distance
      double XDrift = voxel->VoxelID().X() - xyz[0];

      if(XDrift < 0.) continue; 
      
      // Drift time (nano-sec)
      double TDrift = XDrift/driftvelocity;
      double lifetimecorrection = TMath::Exp(-(TDrift/1000.)/electronlifetime);
      double nElectrons = fGeV2Elect * Energy * recomb * lifetimecorrection;
      // Longitudinal & transverse diffusion sigma (cm)
      double LDiffSig = TMath::Sqrt(fLongDiff*TDrift);
      double TDiffSig = TMath::Sqrt(fTranDiff*TDrift);
      int nClus = 1 + (int)(nElectrons/fClusterSize);
      
//       std::cout << "there are " << nClus << " clusters for voxel " << i << std::endl;
      
      // Drift nClus electron clusters to the induction plane
      for(int k = 0; k<nClus; ++k){
	double XDiff = gauss.fire(0.,LDiffSig);
	// Correct drift time for longitudinal diffusion
	double TDiff = TDrift + XDiff/driftvelocity;
	// Smear the Y,Z position by the transverse diffusion
	double YDiff = gauss.fire(voxel->VoxelID().Y(),TDiffSig);
	double ZDiff = gauss.fire(voxel->VoxelID().Z(),TDiffSig);
	
	fDiffuseX->Fill(XDiff);
	fDiffuseY->Fill(voxel->VoxelID().Y() - YDiff);
	fDiffuseZ->Fill(voxel->VoxelID().Z() - ZDiff);

	///grab the nearest channel to the xyz position
	xyz1[0] = xyz[0] + XDiff;
	xyz1[1] = YDiff;
	xyz1[2] = ZDiff;
	
	double nElDiff = 0.;
	if(k<nClus-1) {
	  nElDiff = fClusterSize;
	} 
	else {
	  // Drift the last cluster with smaller size
	  nElDiff = nElectrons - (nClus-1)*fClusterSize;
	}
	
	///make a collection of electrons for each plane
	for(int p = 0; p < planes; ++p){
	  xyz1[0] = xyz[0] + pitch[p];
	  double pitchT = fabs(pitch[p])/driftvelocity;
	  unsigned int channel = geom->NearestChannel(xyz1);
	  fChannels->Fill(channel);
	  double voxXYZ[3] = {voxel->VoxelID().X(), voxel->VoxelID().Y(), voxel->VoxelID().Z()};
	  ecol->push_back(sim::Electrons(channel, TDiff+pitchT, nElDiff, voxel, voxXYZ));
	}//end loop over planes
      }//end loop over clusters
    }//end loop over voxels in list

    ///still save the vector even if there are no electrons
    if(ecol->size() == 0){
      std::cerr << "no electrons made for this event" << std::endl;
    }

    ///sort the electron vector by channel number
    std::sort(ecol->begin(), ecol->end(), chan_sort);

    evt.put(ecol);

    return;

  }


}
