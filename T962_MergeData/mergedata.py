# Configuration file for mergedata
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
#  /author kinga.partyka@yale.edu
#          
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as mergedata

# Give this job a name.  
process = mergedata.Process("MergeData")

# Maximum number of events to do.
#put 8119 for R728
#put 30476 for R729
process.maxEvents = mergedata.untracked.PSet(
    input = mergedata.untracked.int32(16424) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = mergedata.Service(
    "TFileService",
    fileName = mergedata.string("mergedata_hist.root"),
    closeFileFast = mergedata.untracked.bool(False)
)

process.Timing = mergedata.Service("Timing");

# Define the geometry.
process.Geometry = mergedata.Service(
    "Geometry",
    SurfaceY = mergedata.double(130.0e2), # in cm
    Name     = mergedata.string("argoneut"),
    GDML     = mergedata.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Service to get my MC events, which were run up through DetSim.
process.source = mergedata.Source("PoolSource",
                                fileNames = 
                                mergedata.untracked.vstring("/argoneut/data/rootfiles_ART/R609_D20090915_T153438.root")
                                #mergedata.untracked.vstring("/argoneut/app/users/soderber/larsoft_svn/Kinga_Events/testARTevents.root")
                                #skipEvents=mergedata.untracked.uint32(29997) for run 728
                                #remember you need a comma if you use skipevents
                                )


process.merge = mergedata.EDProducer(
    "MergeData",
    daq   = mergedata.string("source") 
   # file       = mergedata.double(1.0),
    
    )    


# Write the events to the output file.
process.output = mergedata.OutputModule(
    "PoolOutputModule",
    fileName = mergedata.untracked.string('mergedata_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = mergedata.EndPath(process.merge*process.output )

