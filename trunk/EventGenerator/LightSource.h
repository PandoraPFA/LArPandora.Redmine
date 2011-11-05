// LightSource.h -  Ben Jones, MIT 2010
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
//  bool    ExtendLibrary       - if library building, extend existing data file or build new library
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



#ifndef EVGEN_LIGHTSOURCE_H
#define EVGEN_LIGHTSOURCE_H
#include <iostream>

#include "art/Framework/Core/EDProducer.h"

#include "fstream"
#include "string"
#include "TLorentzVector.h"
#include "Simulation/PhotonVoxels.h"
#include "TTree.h"

namespace simb { class MCTruth;   }

namespace sim { class PhotonVoxelDef; }

namespace evgen {
  class SingleParticle;

  /// A module for optical MC testing and library building
  class LightSource : public art::EDProducer {
  public:
    explicit LightSource(fhicl::ParameterSet const& pset);
    virtual ~LightSource();                        
  
    void produce(art::Event & evt);
    void beginRun(art::Run& run);

  private:

    void Sample(simb::MCTruth &truth);
         
    int               fSeed;              //random number seed
    std::string       fVersion;           //version of the configuration

    // Flags to mark module modes 
    static const int  kUNIF = 0;
    static const int  kGAUS = 1;
    static const int  kFILE = 0;
    static const int  kSCAN = 1;

    // File stream, filename and empty string for file processing
    ifstream           fInputFile;
    std::string        fFileName;      
    char               fDummyString[256];

    // A ttree to keep track of where particles have been shot - ends up in histos.root
    TTree *            fPhotonsGenerated;
    TLorentzVector     fShotPos;
    TLorentzVector     fShotMom;
    Int_t              fEvID;

    // Parameters loaded from config - both modes
    int                fSourceMode;     // Mode to run in - scan or file
    bool               fFillTree;       // Do we want to create a TTree of shot particles?
    bool               fBuildLibrary;   // Store voxel information in event?
    bool               fExtendLibrary;  // Extend existing library or build new?

    int                fPosDist;        //
    int                fTDist;          // Random distributions to use : 1= gauss, 0= uniform
    int                fPDist;          //  

    //Scan mode specific parameters
    int fXSteps;                        //
    int fYSteps;                        //  Number of steps to take in each dimension
    int fZSteps;                        // 

    sim::PhotonVoxelDef fThePhotonVoxelDef;  // The photon voxel definition object for scan mode

    int fVoxelCount;                    // Total number of voxels
    int fCurrentVoxel;                  // Counter to keep track of vox ID


    //  TPC Measurements
    TVector3 fTPCCenter;    
    TVector3 fTPCDimensions;           
    std::vector<double> fRegionMin;
    std::vector<double> fRegionMax;
    bool fUseCustomRegion;
    

    // Parameters used to shoot in distributions
    double              fX;              // central x position of source
    double              fY;              // central y position of source
    double              fZ;              // central z position of source
    double              fT;              // central t position of source
    double              fSigmaX;         // x width
    double              fSigmaY;         // y width
    double              fSigmaZ;         // z width
    double              fSigmaT;         // t width
    double              fP;              // central momentm of photon 
    double              fSigmaP;         // mom width;

    // Number of photons per event
    int                fN;              // number of photons per event

  };
};

#endif // EVGEN_LIGHTSOURCE_H
////////////////////////////////////////////////////////////////////////
