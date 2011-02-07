#
# $Id: GNUmakefile,v 1.4 2009/06/05 18:49:52 t962cvs Exp $
#
# messier@indiana.edu
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := VertexActivity
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lUtilities -lGeometry -lSimulationBase -lSimulation -lRawData -lRecoBase  -lFilters  -l$(PACKAGE) 


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk

override LIBLIBS += -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service
