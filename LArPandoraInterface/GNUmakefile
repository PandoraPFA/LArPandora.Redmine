#
# Makefile for LArPandoraInterface
#
include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

#
#ROOTLIBS  = $(shell root-config --libs) -lEG -lTreePlayer -lGeom -lFFTW -lReflex -lEGPythia6 -lGui -lEve
#
#override CPPFLAGS += $(shell root-config --cflags)
#override LDFLAGS += $(shell root-config --libs) 
#
#override LOADLIBES += $(ROOTLIBS) $(ROOTGLIBS) -lRGL -lGed -lX3d
#

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(NUTOOLS_LIB) -lSimulationBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lGeometry -lSimulation -lRawData -lRecoBase -lUtilities -lLArPandoraAlgorithms -L${PANDORA_LIB} -lPandoraSDK

override CPPFLAGS += -I${NUTOOLS_INC} -I${PANDORA_INC}


