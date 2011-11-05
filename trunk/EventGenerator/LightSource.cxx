// LightSource.cxx  - Ben Jones, MIT 2010
//
// Light source event generator which simulate an extended isotropic photon source
//
// The light source can be run in two modes, file mode or scan mode.  Each requires
// the specification of a different set of parameters.
//
// FILE MODE :
// Light source position, intensity and shape are supplied on an event by event basis
// in a text file.  See the example provided for the format. Pararmeters required:
//
//  int32   SourceMode = 0      - sets light source to file mode
//  string  FileName            - file of per event light source specifications
//  int32   PosDist             - how to distribute production points sampled in momentum, position
//  int32   PDist                   and time ranges specified.  For all of these :
//  int32   TDist                   0 = uniform and 1 = gauss
//  bool    FillTree            - whether to write a tree of photon production points to fileservice
//  bool    BuildLibrary        - whether to write the light source specification to the event for library building / analysis
//  bool    ExtendLibrary       - if building a library, extend existing data file or build a new one
//
// Upon reaching the end of the file, the light source will loop back to the first point.
// hence a one line text file will give a constant light source size, position and intensity.
//
// SCAN MODE:
// Divide volume into cuboidal regions and produce an isotropic light source in each,
// using one region per event.  User can specify either to use the full detector volume
// or some custom specified volume.
//
// This mode is used when building a fast photon sim library, and performing volume
// scan sensitivity studies.
//
//  int32   SourceMode = 1      - sets light source to scan mode
//  int32   N                   - number of photons to shoot from each point
//  double  P                   - peak photon momentum (or energy) in eV
//  double  SigmaP              - momentum distribution width
//  double  XSteps              - Number of regions to divide volume into in each direction 
//  double  YSteps
//  double  ZSteps
//  double  T0                  - Peak time of photon production
//  double  SigmaT              - time distribution width
//  int32   PosDist             - how to distribute production points sampled in momentum, position
//  int32   PDist                 and time ranges specified.  For all of these :
//  int32   TDist                   0 = uniform and 1 = gaussian
//  bool    FillTree            - whether to write a tree of photon production points to fileservice
//  bool    BuildLibrary        - whether to write the light source specificationto the event for library building / analysis
//  bool    UseCustomRegion     - supply our own volme specification or use the full detector volume?
//  vdouble[3]  RegionMin       - bounding corners of the custom volume specification 
//  vdouble[3]  RegionMax           (only used if UseCustomRegion=true)
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <memory>


// ART includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools includes
#include "SimulationBase/simbase.h"
#include "EventGeneratorBase/evgenbase.h"

// lar includes
#include "EventGenerator/LightSource.h"
#include "Geometry/geo.h"
#include "Simulation/PhotonVoxels.h"
#include "Simulation/PhotonLibraryParameters.h"
#include "PhotonPropagation/PhotonLibraryService.h"
#include "SummaryData/summary.h"

#include "TVector3.h"
#include "TDatabasePDG.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace evgen{
  
  //----------------------------------------------------------------
  LightSource::LightSource(fhicl::ParameterSet const& pset) 
  {

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
 
    fSourceMode   =     (pset.get<int >("SourceMode")  );
    fFillTree     =     (pset.get<bool>("FillTree")    );
    fBuildLibrary =     (pset.get<bool>("BuildLibrary"));
    fPosDist      =     (pset.get<int >("PosDist")     );
    fPDist        =     (pset.get<int >("PDist")       );
    fTDist        =     (pset.get<int >("TDist")       );
   unsigned int seed = pset.get< unsigned int >("Seed", evgb::GetRandomNumberSeed());

    createEngine(seed);

    // load optional parameters in function
    produces< sumdata::RunData, art::InRun >();
    produces< std::vector<simb::MCTruth> >();
    if(fBuildLibrary) produces<std::vector<sim::PhotonLibraryParameters> >();

    if(fSourceMode==kFILE)
      {
	fFileName  = pset.get<std::string>("SteeringFile");
	fInputFile.open(fFileName.c_str());
	fInputFile.getline(fDummyString,256);

      }
    else if (fSourceMode==kSCAN)
      {
	fXSteps = pset.get<int >("XSteps");
	fYSteps = pset.get<int >("YSteps");
	fZSteps = pset.get<int >("ZSteps");

	fT      = pset.get<double>("T0");
	fSigmaT = pset.get<double>("SigmaT");
	fN      = pset.get<int   >("N");
	
	fP      = pset.get<double>("P");
	fSigmaP = pset.get<double>("SigmaP");

	fUseCustomRegion = pset.get<bool>("UseCustomRegion");

	if(fUseCustomRegion)
	  {	
	    fRegionMin = pset.get< std::vector<double> >("RegionMin");
	    fRegionMax = pset.get< std::vector<double> >("RegionMax");
	  }

	// Get dimensions of TPC volume
	art::ServiceHandle<geo::Geometry> geo;
	fTPCCenter     = geo->GetTPCFrontFaceCenter();
	fTPCDimensions = TVector3(geo->DetHalfWidth(),
				  geo->DetHalfHeight(),
				  geo->DetLength());
	
	std::cout<<"Light Source Reading TPC Dimensions" << std::endl<<
	  "  Half Width " << fTPCDimensions[0] << "   Half Height " <<fTPCDimensions[1] << "   Length " << fTPCDimensions[2] <<std::endl <<
	  "  Centre X " << fTPCCenter[0] << "  Center Y " << fTPCCenter[1] << "  CenterZ " << fTPCCenter[2] <<std::endl;
	
	

	fCurrentVoxel=0;

	std::cout<<"Light Source Debug : fBuildLibrary " << fBuildLibrary<<std::endl;


	// define voxelization based on parameters read from config.
	// There are two modes - either read the dimensions of the TPC from
	// the geometry, or use values specified by the user.
	if(!fUseCustomRegion)
	fThePhotonVoxelDef = sim::PhotonVoxelDef(fTPCCenter[0]-fTPCDimensions[0],     // min X
						fTPCCenter[0]+fTPCDimensions[0],     // max x
						fXSteps,                             // steps in x
						fTPCCenter[1]-fTPCDimensions[1],     // min y
						fTPCCenter[1]+fTPCDimensions[1],     // max y
						fYSteps,                             // steps in y
						fTPCCenter[2],                       // min z
						fTPCCenter[2]+fTPCDimensions[2],     // max z
						fZSteps                              // steps in z
						);
	else
	  {
	    fThePhotonVoxelDef = sim::PhotonVoxelDef(fRegionMin[0], 
							fRegionMax[0],
							fXSteps,
							fRegionMin[1],
							fRegionMax[1],
							fYSteps,
							fRegionMin[2],
							fRegionMax[2],
							fZSteps);   
	  }

	// If this a library building job, initiate the library service
	if(fBuildLibrary)
	  {
	    fExtendLibrary = pset.get<bool>("ExtendLibrary");
	    art::ServiceHandle<phot::PhotonLibraryService> pls;
	    pls->InitializeBuildJob(fThePhotonVoxelDef, fExtendLibrary);
	  }

	// Set distribution widths to voxel size

	fSigmaX = fThePhotonVoxelDef.GetVoxelSize().X()/2.0;
	fSigmaY = fThePhotonVoxelDef.GetVoxelSize().Y()/2.0;
	fSigmaZ = fThePhotonVoxelDef.GetVoxelSize().Z()/2.0;

	// Get number of voxels we will step through
	
	fVoxelCount = fThePhotonVoxelDef.GetNVoxels();

	std::cout<<"Light Source : Determining voxel params : " <<fVoxelCount<<" "<< fSigmaX<<" "<< fSigmaY<<" " <<fSigmaZ<<std::endl; 

      }
    else
      {
	std::cout<<"EVGEN Light Source : Unrecognised light source mode" <<std::endl;
	assert(0);
      }
    



    if(fFillTree)
      {
	art::ServiceHandle<art::TFileService> tfs;
	fPhotonsGenerated = tfs->make<TTree>("PhotonsGenerated","PhotonsGenerated");
	fPhotonsGenerated->Branch("X",&(fShotPos[0]),"X/D");
	fPhotonsGenerated->Branch("Y",&(fShotPos[1]),"Y/D");
	fPhotonsGenerated->Branch("Z",&(fShotPos[2]),"Z/D");
	fPhotonsGenerated->Branch("T",&(fShotPos[3]),"T/D");
	fPhotonsGenerated->Branch("PX",&(fShotMom[0]),"PX/D");
	fPhotonsGenerated->Branch("PY",&(fShotMom[1]),"PY/D");
	fPhotonsGenerated->Branch("PZ",&(fShotMom[2]),"PZ/D");
	fPhotonsGenerated->Branch("PT",&(fShotMom[3]),"PT/D");
	fPhotonsGenerated->Branch("EventID",&fEvID,"EventID/I");
      }
  }


  //----------------------------------------------------------------
  LightSource::~LightSource()
  {
  }

  //____________________________________________________________________________
  void LightSource::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;

    geo::DetId_t detid = geo->DetId();

    std::auto_ptr<sumdata::RunData> runcol(new sumdata::RunData(detid));

    run.put(runcol);

    return;
  }

  //----------------------------------------------------------------
  void LightSource::produce(art::Event& evt)
  {
    
    // FILE MODE -
    //  Each event, read coordinates of gun and number of photons to shoot from file
    
    if(fSourceMode==kFILE)
      {
	// Loop file if required
	if(fInputFile.eof())
	  {
	    std::cout<<"EVGEN Light Source : Warning, reached end of file, looping back to beginning"<<std::endl;
	    fInputFile.seekg(0,std::ios::beg);
	    fInputFile.clear();
	  }
	
	if(!fInputFile.is_open() | fInputFile.fail() )
	  {
	    std::cout<<"EVGEN Light Source : File error in " <<fFileName.c_str()<<std::endl;
	    assert(0);
	  }
	else
	  { 
	    // read in one line
	    fInputFile >> fX >> fY >> fZ >> fT >> fSigmaX >> fSigmaY >> fSigmaZ >> fSigmaT >> fP >> fSigmaP >> fN;
	    fInputFile.getline(fDummyString,256);
	    fThePhotonVoxelDef = sim::PhotonVoxelDef(fX - fSigmaX, 
							fX + fSigmaX,
							1,
							fY - fSigmaY,
							fY + fSigmaY,
							1,
							fZ - fSigmaZ,
							fZ + fSigmaZ,
						        1);   
	    
	    fCurrentVoxel=0;
	  }
      }


    // SCAN MODE -
    //  Step through detector using a number of steps provided in the config file
    //  firing a constant number of photons from each point
    else if(fSourceMode==kSCAN)
      {
	if(fCurrentVoxel!=fVoxelCount) 
	  {
	    fCurrentVoxel++;
	  }
	else
	  {
	    std::cout<<"EVGEN Light Source fully scanned detector.  Starting over."<<std::endl;
	    fCurrentVoxel=0;
	  }

	TVector3 VoxelCenter = fThePhotonVoxelDef.GetPhotonVoxel(fCurrentVoxel).GetCenter();
	
	fX = VoxelCenter.X();
	fY = VoxelCenter.Y();
	fZ = VoxelCenter.Z();
	
      }
	
    
    // UNRECOGNISED MODE 
    //  - neither file or scan mode, probably a config file error

    else
      {
	std::cout<<"EVGEN : Light Source, unrecognised source mode" <<std::endl;
	assert(0);
      }
    


    //    std::cout<<"EVGEN Light source to be placed at  (x, y, z, t, dx, dy, dz, dt, n) " <<
    //      fX << " " << fY << " " <<fZ << " " << fT << " " << 
    //      fSigmaX << " " << fSigmaY << " " << fSigmaZ << " " << fSigmaT << " " <<
    //      fN << std::endl<<std::endl;
    
      
    std::auto_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    Sample(truth);
    
    //     std::cout << "put mctruth into the vector" << std::endl;
    truthcol->push_back(truth);

    //     std::cout << "add vector to the event " << truthcol->size() << std::endl;
    evt.put(truthcol);

    if(fBuildLibrary)
      {
	std::cout<<"Light source : Storing voxel params in event"<<std::endl;
	sim::PhotonLibraryParameters LibraryParams(fThePhotonVoxelDef, fCurrentVoxel,fN);
	std::auto_ptr< std::vector<sim::PhotonLibraryParameters> > paramcol(new std::vector<sim::PhotonLibraryParameters>);
	paramcol->push_back(LibraryParams);
	evt.put(paramcol);
      }

    return;


  }


  void LightSource::Sample(simb::MCTruth& mct) 
  {
    std::cout<<"Light source debug : Shooting at " << fX <<" " << fY<<" "<< fZ<<std::endl;
    
    // get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat   flat(engine);
    CLHEP::RandGaussQ gauss(engine);

    for(int j=0; j!=fN; j++)
      {
	// Choose momentum (supplied in eV, convert to GeV)
	double p = fP;
	if (fPDist == kGAUS) {
	  p = gauss.fire(fP, fSigmaP);
	}
	else {
	  p = fP + fSigmaP*(2.0*flat.fire()-1.0);
	}
	p /= 1000000000.;
	
	// Choose position
	TVector3 x;
	if (fPosDist == kGAUS) {
	  x[0] = gauss.fire(fX, fSigmaX);
	  x[1] = gauss.fire(fY, fSigmaY);
	  x[2] = gauss.fire(fZ, fSigmaZ);
	}
	else {
	  x[0] = fX + fSigmaX*(2.0*flat.fire()-1.0);
	  x[1] = fY + fSigmaY*(2.0*flat.fire()-1.0);
	  x[2] = fZ + fSigmaZ*(2.0*flat.fire()-1.0);
	}
	
	// Choose time
	double t;
	if (fTDist == kGAUS) {
	  t = gauss.fire(fT, fSigmaT);
	}
	else {
	  t = fT + fSigmaT * (2.0 * flat.fire()-1.0);
	}
	
	
	//assume the position is relative to the center of the TPC
	x += fTPCCenter;
	
	fShotPos = TLorentzVector(x[0], x[1], x[2], t);
	
	//	std::cout << "Light source placing photon (x,y,z,t) " << x[0] << " " << x[1] << " " << x[2] << " " << t <<std::endl;
	
	
	double thxz,thyz;

	// Choose angles
	thxz = 180 * flat.fire();
        thyz = 360 * flat.fire();
	
	// Generate momentum 4-vector
	fShotMom = TLorentzVector(0.0,0.0,p,p);
	fShotMom.RotateY(thxz*M_PI/180.0);
	fShotMom.RotateX(thyz*M_PI/180.0);
	
	//	std::cout<<"Light source chose particle direction (th_xz, th_yz) " << thxz  << ", " << thyz <<std::endl;
      

       
	int trackid = -1*(j+1); // set track id to -i as these are all primary particles and have id <= 0
	std::string primary("primary");
	int PDG=0; //optical photons have PDG 0

	simb::MCParticle part(trackid, PDG, primary);
	part.AddTrajectoryPoint(fShotPos, fShotMom);

	if(fFillTree)
	  fPhotonsGenerated->Fill();
	
	mct.Add(part);
      }
  }

}
