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

override LOADLIBES += $(ROOTLIBS)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lGeometry -lSimulationBase -lSimulation 
