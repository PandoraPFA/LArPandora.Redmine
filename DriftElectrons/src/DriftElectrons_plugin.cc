////////////////////////////////////////////////////////////////////////
/// \file  DriftElectrons.cxx
/// \brief Module to drift ionization electrons to wire planes
///
/// \version $Id: DriftElectrons.cxx,v 1.20 2010/04/30 16:10:37 brebel Exp $
/// \author  baller@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DriftElectrons/inc/DriftElectrons.h"
#include "Geometry/inc/geo.h"
#include "Simulation/inc/LArVoxelList.h"
#include "Simulation/inc/LArVoxelData.h"
#include "Simulation/inc/Electrons.h"

#include "TH1.h"
#include "TMath.h"
#include "TRandom3.h"
#include <ctime>

namespace dfe{

  // A macro required for a JobControl module.
  DEFINE_FWK_MODULE(DriftElectrons);

  static bool chan_sort(const sim::Electrons e1, const sim::Electrons e2)
  {
    return e1.Channel() < e2.Channel();
  }

  //-----------------------------------------------------
  DriftElectrons::DriftElectrons(edm::ParameterSet const& pset) :
    fRecombFactor    (pset.getParameter< double >("RecombinationFactor")),
    fLongDiff        (pset.getParameter< double >("LongDiff")),
    fTranDiff        (pset.getParameter< double >("TranDiff")),
    fDriftVel        (pset.getParameter< double >("DriftVel")),
    fClusterSize     (pset.getParameter< double >("ClusterSize")),
    fLArG4ModuleLabel(pset.getParameter< std::string >("LArG4ModuleLabel"))
  {
    produces< std::vector<sim::Electrons> >();

    // Conversion from GeV energy loss to electrons
    // 0.7 = recombination factor for 500 V/cm
    // 23.6 eV per ion pair * 1E9 ev/GeV
    fGeV2Elect = fRecombFactor * 1.E9/23.6 ;

    // There's a factor of two in the diffusion equation
    fLongDiff *= 2.;
    fTranDiff *= 2.;

    int RanSeed = time(0);
    fRandom.SetSeed(RanSeed);
  }

  //-----------------------------------------------------
  DriftElectrons::~DriftElectrons()
  {
  }

  //-----------------------------------------------------
  void DriftElectrons::beginJob(edm::EventSetup const&)
  {
    // get access to the TFile service
    edm::Service<edm::TFileService> tfs;

    edm::Service<geo::Geometry> geo;

    fChannels    = tfs->make<TH1D>("channels",  ";channel;# electrons",  geo->Nchannels(), 0, geo->Nchannels());
    fDiffuseX    = tfs->make<TH1D>("diffusex",  ";x diffusion (cm);",    200, -5., 5.  );
    fDiffuseY    = tfs->make<TH1D>("diffusey",  ";y diffusion (cm);",    200, -5., 5.  );
    fDiffuseZ    = tfs->make<TH1D>("diffusez",  ";z diffusion (cm);",    200, -5., 5.  );
    fNumVoxels   = tfs->make<TH1D>("numvoxels", ";voxels in event;",     500,  0., 5000.);
    fVoxelEnergy = tfs->make<TH1D>("voxele",    ";energy (GeV); voxels", 500,  0., 0.01  );
  }

  //-----------------------------------------------------
  void DriftElectrons::produce(edm::Event& evt, edm::EventSetup const&)
  {
    
    //get a collection of electrons
    std::auto_ptr<std::vector<sim::Electrons> > ecol(new std::vector<sim::Electrons>);

    edm::Service<geo::Geometry> geom;
  
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
    ///reset the xyz array to be at the first induction plane
    geom->Plane(0).Wire(0).GetCenter(xyz);

//     std::cout << "In Drift: X =" << xyz[0] << std::endl;
  
    // Read in the LArVoxelList object(s).
    edm::Handle< std::vector<sim::LArVoxelData> > vxlistHandle;
    evt.getByLabel(fLArG4ModuleLabel,vxlistHandle);
    const std::vector<sim::LArVoxelData> &larVoxelList(*vxlistHandle);

    fNumVoxels->Fill(larVoxelList.size());

    // There's probably only one LArVoxelList per event, but FMWK
    // always reads a vector of pointers.  For each LArVoxelList:
    for(unsigned int i = 0; i < larVoxelList.size(); ++i){

      // Get the reference to the LArVoxelID in the LArVoxelList.
      edm::Ptr<sim::LArVoxelData> voxel(vxlistHandle, i);
      
      // Get the energy calculated by the LArVoxelList method.
      double Energy = voxel->Energy();
      fVoxelEnergy->Fill(Energy);
      
      double nElectrons = fGeV2Elect * Energy;
      // X drift distance
      double XDrift = voxel->VoxelID().X() - xyz[0];
      
//       std::cout << XDrift << " " << xyz[0] << " " << voxel->VoxelID().X() << std::endl;
      if(XDrift < 0.) continue; 
      // Drift time (nano-sec)
      double TDrift = XDrift/fDriftVel;
      // Longitudinal & transverse diffusion sigma (cm)
      double LDiffSig = TMath::Sqrt(fLongDiff*TDrift);
      double TDiffSig = TMath::Sqrt(fTranDiff*TDrift);
      int nClus = 1 + (int)(nElectrons/fClusterSize);
      
//       std::cout << "there are " << nClus << " clusters for voxel " << i << std::endl;
      
      // Drift nClus electron clusters to the induction plane
      for(int k = 0; k<nClus; ++k){
	double XDiff = fRandom.Gaus(0.,LDiffSig);
	// Correct drift time for longitudinal diffusion
	double TDiff = TDrift + XDiff/fDriftVel;
	// Smear the Y,Z position by the transverse diffusion
	double YDiff = fRandom.Gaus(voxel->VoxelID().Y(),TDiffSig);
	double ZDiff = fRandom.Gaus(voxel->VoxelID().Z(),TDiffSig);
	
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
	  double pitchT = fabs(pitch[p])/fDriftVel;
	  unsigned int channel = geom->NearestChannel(xyz1);
	  fChannels->Fill(channel);
	  ecol->push_back(sim::Electrons(channel, TDiff+pitchT, nElDiff, voxel));
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
