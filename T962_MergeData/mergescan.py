# Configuration file for mergescan
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
#  /author joshua.spitz@yale.edu
#          
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mergescan

# Give this job a name.  
process = mergescan.Process("MergeScan")

# Maximum number of events to do.
process.maxEvents = mergescan.untracked.PSet(
    input = mergescan.untracked.int32(3) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mergescan.Service(
    "TFileService",
    fileName = mergescan.string("mergescan_hist.root"),
    closeFileFast = mergescan.untracked.bool(False)
)

process.Timing = mergescan.Service("Timing");

# Define the geometry.
process.Geometry = mergescan.Service(
    "Geometry",
    SurfaceY = mergescan.double(130.0e2), # in cm
    Name     = mergescan.string("argoneut"),
    GDML     = mergescan.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Service to get my MC events, which were run up through DetSim.
process.source = mergescan.Source("PoolSource",
                                fileNames = mergescan.untracked.vstring("/argoneut/app/users/soderber/larsoft_svn/Kinga_Events/testARTevents.root")
                                )


process.merge = mergescan.EDProducer(
    "MergeScan",
    daq   = mergescan.string("source"),
    scanners = mergescan.vstring("spitz7","patch")
    )    


# Write the events to the output file.
process.output = mergescan.OutputModule(
    "PoolOutputModule",
    fileName = mergescan.untracked.string('mergescan_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = mergescan.EndPath(process.merge*process.output )

