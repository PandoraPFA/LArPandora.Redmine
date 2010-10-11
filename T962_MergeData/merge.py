# Configuration file for drifting electrons
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author B. Rebel
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as mergedata

# Give this job a name.  
process = mergedata.Process("MergeData")

# Maximum number of events to do.
process.maxEvents = mergedata.untracked.PSet(
    input = mergedata.untracked.int32(10)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mergedata.Service("TFileService",
                       fileName = mergedata.string("mergedata_hist.root"),
                       closeFileFast = mergedata.untracked.bool(False)
)


# Initialize the random number sequences.
#process.add_(mu2e.Service("RandomNumberService",
#                          globalSeed=mu2e.untracked.int32(9877)
#))

# Define the geometry.
process.Geometry = mergedata.Service(
    "Geometry",
    SurfaceY = mergedata.double(130.0e2), # in cm
    Name     = mergedata.string("argoneut"),
    GDML     = mergedata.FileInPath("Geometry/gdml/argoneut.gdml")
)

# Define the service for passing voxel parameters
process.MergeData = mergedata.Service(
    "MergeData",
    VoxelSizeX               = mergedata.double(0.03),   #in cm
    VoxelSizeY               = mergedata.double(0.03),   #in cm
    VoxelSizeZ               = mergedata.double(0.03),   #in cm
    VoxelSizeT               = mergedata.double(5000.0), #in ns
    VoxelOffsetX             = mergedata.double(0.0),    #in cm
    VoxelOffsetY             = mergedata.double(0.0),    #in cm
    VoxelOffsetZ             = mergedata.double(0.0),    #in cm
    VoxelOffsetT             = mergedata.double(-2500.0),#in ns
    VoxelEnergyCut           = mergedata.double(1.e-6)   #in GeV
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = mergedata.Source(
    "PoolSource",
    fileNames = mergedata.untracked.vstring("larg4.root")
)

# Make some particles to go through
process.mergedatal = mergedata.EDProducer(
    "mergedatalectrons",
    LArG4ModuleLabel         = mergedata.string("largeant"),
    LongDiff                 = mergedata.double(6.2e-9), #in cm^2/ns
    TranDiff                 = mergedata.double(16.3e-9),#in cm^2/ns
    DriftVel                 = mergedata.double(1.55e-4),#in cm/ns
    RecombinationFactor      = mergedata.double(0.7),
    ClusterSize              = mergedata.double(600.0)
)

# Write the events to the output file.
process.output = mergedata.OutputModule(
    "PoolOutputModule",
    fileName = mergedata.untracked.string('file:mergedata.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths.
process.doit = mergedata.EndPath( process.mergedatal*process.output )

