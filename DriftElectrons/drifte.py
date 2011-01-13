# Configuration file for drifting electrons
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author B. Rebel
#
# Spacing is not signficant in this file.

# Define the default configuratio for the framework.
import FWCore.ParameterSet.python.Config as drifte

# Give this job a name.  
process = drifte.Process("DriftElectrons")

# Maximum number of events to do.
process.maxEvents = drifte.untracked.PSet(
    input = drifte.untracked.int32(10)
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = drifte.Service("TFileService",
                       fileName = drifte.string("drifte_hist.root"),
                       closeFileFast = drifte.untracked.bool(False)
)


# Initialize the random number sequences.
#process.add_(mu2e.Service("RandomNumberService",
#                          globalSeed=mu2e.untracked.int32(9877)
#))

# Define the geometry.
process.Geometry = drifte.Service(
    "Geometry",
    SurfaceY = drifte.double(130.0e2), # in cm
    Name     = drifte.string("argoneut"),
    GDML     = drifte.FileInPath("Geometry/gdml/argoneut.gdml")
)

# Define the LAr property service
process.LArProperties = drifte.Service(
    "LArProperties",
    Efield       = drifte.double(0.5),   #kV/cm
    Temperature  = drifte.double(87.),     #kelvin
    Electronlifetime   = drifte.double(1000000000.)     #microseconds
)

# Define the service for passing voxel parameters
process.LArVoxelCalculator = drifte.Service(
    "LArVoxelCalculator",
    VoxelSizeX               = drifte.double(0.03),   #in cm
    VoxelSizeY               = drifte.double(0.03),   #in cm
    VoxelSizeZ               = drifte.double(0.03),   #in cm
    VoxelSizeT               = drifte.double(5000.0), #in ns
    VoxelOffsetX             = drifte.double(0.0),    #in cm
    VoxelOffsetY             = drifte.double(0.0),    #in cm
    VoxelOffsetZ             = drifte.double(0.0),    #in cm
    VoxelOffsetT             = drifte.double(-2500.0),#in ns
    VoxelEnergyCut           = drifte.double(1.e-6)   #in GeV
)

# Define and configure some modules to do work on each event.
# Modules are just defined for now, the are scheduled later.

# Start each new event with an empty event.
process.source = drifte.Source(
    "PoolSource",
    fileNames = drifte.untracked.vstring("larg4.root")
)

# Make some particles to go through
process.driftel = drifte.EDProducer(
    "DriftElectrons",
    LArG4ModuleLabel         = drifte.string("largeant"),
    LongDiff                 = drifte.double(6.2e-9), #in cm^2/ns
    TranDiff                 = drifte.double(16.3e-9),#in cm^2/ns
    RecombA                  = drifte.double(0.8),
    Recombk                  = drifte.double(0.097),  #in g/(MeVcm^{2})
    ClusterSize              = drifte.double(600.0)
)

# Write the events to the output file.
process.output = drifte.OutputModule(
    "PoolOutputModule",
    fileName = drifte.untracked.string('file:drifte.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths.
process.doit = drifte.EndPath( process.driftel*process.output )

