#
# $Id: GNUmakefile,v 1.8 2009/11/20 19:04:17 soderber Exp $
#
# brebel@fnal.gov
#
#
ROOTCINT    := YES
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

ROOTLIBS  = $(shell root-config --libs) -lEG -lTreePlayer -lGeom -lFFTW -lReflex -lEGPythia6 -lGui

override CPPFLAGS += $(shell root-config --cflags) 
override LDFLAGS += $(shell root-config --libs) 

override LOADLIBES += $(ROOTLIBS) $(ROOTGLIBS) -lRGL -lGed -lX3d

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk

override LIBLIBS += $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lFilters -lGeometry -lSimulationBase -lSimulation -lRawData -lRecoBase -lEventDisplayBase -lEventDisplay -lUtilities -lT962_Objects

