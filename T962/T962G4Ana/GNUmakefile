#
# $Id: GNUmakefile,v 1.4 2009/06/19 19:13:26 brebel Exp $
#
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_geant4.mk
include SoftRelTools/arch_spec_xercesc.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

CPPFLAGS     += -I$(G4INSTALL)/include

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
########################################################################

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(CLHEP_BASE)/lib -lCLHEP -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lEventGeneratorBase -lSimulationBase -lSimulation -lGeometry 