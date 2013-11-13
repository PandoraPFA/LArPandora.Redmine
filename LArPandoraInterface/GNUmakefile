#
# Makefile for LArPandoraInterface
#
include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_nutools.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(NUTOOLS_LIB) -lSimulationBase -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lGeometry -lSimulation -lRawData -lRecoBase -lUtilities -lLArPandoraAlgorithms -L$(PANDORA_LIB) -lPandoraSDK -lPandoraMonitoring

override CPPFLAGS += -I$(NUTOOLS_INC) -I$(PANDORA_INC) -I$(PANDORA_FQ_DIR)/PandoraSDK/include -I$(PANDORA_FQ_DIR)/PandoraMonitoring/include -I$(SRT_PRIVATE_CONTEXT)/LArPandoraAlgorithms -I$(SRT_PUBLIC_CONTEXT)/LArPandoraAlgorithms


