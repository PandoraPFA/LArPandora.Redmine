# Configuration file for HarrisVertexFinder
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author J. Spitz   
#   joshua.spitz@yale.edu
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as vertexactivity

# Give this job a name.  
process = vertexactivity.Process("VertexFinder")

# Maximum number of events to do.
process.maxEvents = vertexactivity.untracked.PSet(
    input = vertexactivity.untracked.int32(99) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = vertexactivity.Service(
    "TFileService",
    fileName = vertexactivity.string("vertexactivity_hist.root"),
    closeFileFast = vertexactivity.untracked.bool(False)
)

process.Timing = vertexactivity.Service("Timing");

# Define the geometry.
process.Geometry = vertexactivity.Service(
    "Geometry",
    SurfaceY = vertexactivity.double(130.0e2), # in cm
    Name     = vertexactivity.string("argoneut"),
    GDML     = vertexactivity.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Define the FFT service
process.LArFFT = vertexactivity.Service(
    "LArFFT",
    FFTSize   = vertexactivity.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = vertexactivity.string("ES")
)

# Define the LAr property service
process.LArProperties = vertexactivity.Service(
    "LArProperties",
    Efield       = vertexactivity.double(0.5),   #kV/cm
    Temperature  = vertexactivity.double(87.),     #kelvin
    Electronlifetime   = vertexactivity.double(750.)     #microseconds
)

# Define the service for passing voxel parameters
process.LArVoxelCalculator = vertexactivity.Service(
    "LArVoxelCalculator",
    VoxelSizeX               = vertexactivity.double(0.03),   #in cm
    VoxelSizeY               = vertexactivity.double(0.03),   #in cm
    VoxelSizeZ               = vertexactivity.double(0.03),   #in cm
    VoxelSizeT               = vertexactivity.double(5000.0), #in ns
    VoxelOffsetX             = vertexactivity.double(0.0),    #in cm
    VoxelOffsetY             = vertexactivity.double(0.0),    #in cm
    VoxelOffsetZ             = vertexactivity.double(0.0),    #in cm
    VoxelOffsetT             = vertexactivity.double(-2500.0),#in ns
    VoxelEnergyCut           = vertexactivity.double(1.e-6)   #in GeV
)

# Service to get my MC events, which were run up through DetSim.
process.source = vertexactivity.Source("PoolSource",
                                fileNames = vertexactivity.untracked.vstring("/argoneut/app/users/spitz7/larsoft_9/single_gen.root")
                                )

# process.scanfilter = vertexactivity.EDProducer(
#     "ScanFilter",
#     ScanModuleLabel   = vertexactivity.string("merge"),
#     DigitModuleLabel   = vertexactivity.string("wiresim"),
#     Neutrino_req      = vertexactivity.int32(2),  #0=no neutrinos, 1=maybe neutrino, 2=neutrino. 2 is most stringent, 0 is least stringent  
#     NumShowers_req    = vertexactivity.int32(0),  #Maximum number of showers required to pass
#     NumTracks_req     = vertexactivity.int32(4)   #Maximum number of tracks in any plane required to pass  
#     )


process.caldataCal = vertexactivity.EDProducer(
    "CalWire",
    DigitModuleLabel  = vertexactivity.string("wiresim"),
    ResponseFile       = vertexactivity.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
    ExponentialEndBins = vertexactivity.int32(300),
    UseRawData         = vertexactivity.int32(0)
    )

process.ffthit = vertexactivity.EDProducer(
    "FFTHitFinder",
    CalDataModuleLabel   = vertexactivity.string("caldataCal"),
    MinSigInd       = vertexactivity.double(6.0),
    MinSigCol       = vertexactivity.double(11.0),
    IndWidth        = vertexactivity.double(5.0),
    ColWidth        = vertexactivity.double(7.5),
    Drift           = vertexactivity.double(0.03069),
    # TPC's drift Velocity in units cm/(ADC Sample Time)
    POffset         = vertexactivity.double(20.0),
    OOffset         = vertexactivity.double(24.0),
    MaxMultiHit     = vertexactivity.int32(3)
    )

process.dbscan = vertexactivity.EDProducer(
    "DBcluster",
    HitsModuleLabel   = vertexactivity.string("ffthit"),
    eps               = vertexactivity.double(1.0),
    eps2              = vertexactivity.double(0.75),
    minPts            = vertexactivity.int32(2)
    ) 

process.hough = vertexactivity.EDProducer(
    "HoughLineFinder",
    DBScanModuleLabel        = vertexactivity.string("dbscan"),
    MaxLines                 = vertexactivity.int32(5),
    MinHits                  = vertexactivity.int32(3),
    SaveAccumulator          = vertexactivity.int32(0),
    NumAngleCells            = vertexactivity.int32(10000),
    RhoResolutionFactor      = vertexactivity.int32(10),
    SmootherSigma            = vertexactivity.double(0.),
    MaxDistance              = vertexactivity.double(5.),
    RhoZeroOutRange          = vertexactivity.int32(0),
    ThetaZeroOutRange        = vertexactivity.int32(0),
    HitsPerCluster           = vertexactivity.int32(1)
    )
    
process.harris = vertexactivity.EDProducer(
    "HarrisVertexFinder",
    DBScanModuleLabel    = vertexactivity.string("dbscan"),
    TimeBins             = vertexactivity.int32(256),
    MaxCorners           = vertexactivity.int32(25),
    Gsigma               = vertexactivity.double(1.),
    Window               = vertexactivity.int32(5),
    Threshold            = vertexactivity.double(0.01),
    SaveVertexMap        = vertexactivity.int32(-1)
    )
    
process.vertexmatch = vertexactivity.EDProducer(
    "VertexMatch",    
    HoughModuleLabel       = vertexactivity.string("hough"),
    VertexModuleLabel      = vertexactivity.string("harris"),
    MaxDistance            = vertexactivity.double(30.)
    )  
  
process.vertexactivity = vertexactivity.EDProducer(
    "VertexActivity",    
    DBScanModuleLabel           = vertexactivity.string("dbscan"),
    LArG4ModuleLabel            = vertexactivity.string("largeant"),
    VertexModuleLabel           = vertexactivity.string("vertexmatch"),
    GenieGenModuleLabel         = vertexactivity.string("singlegen"),#geniegen for neutrinos
    ScanModuleLabel             = vertexactivity.string("merge"),
    Cathodetimelocation         = vertexactivity.double(75.),
    Delta_Cathodetimelocation   = vertexactivity.double(10.),
    E_lifetime                  = vertexactivity.double(750.),
    Delta_E_lifetime            = vertexactivity.double(20.),
    Recombination_factor        = vertexactivity.double(0.640),
    Delta_Recombination_factor  = vertexactivity.double(0.013),
    Workfunction_factor         = vertexactivity.double(23.6),
    Delta_Workfunction_factor   = vertexactivity.double(0.5),
    Calibration_factor          = vertexactivity.double(.133),#1/7.54
    Delta_Calibration_factor    = vertexactivity.double(0.),
    ActivityRadius              = vertexactivity.double(25.)#in units of wires
    )     

# Write the events to the output file.
process.output = vertexactivity.OutputModule(
    "PoolOutputModule",
    fileName = vertexactivity.untracked.string('vertexactivity_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = vertexactivity.EndPath(process.caldataCal*process.ffthit*process.dbscan*process.vertexactivity*process.output )

