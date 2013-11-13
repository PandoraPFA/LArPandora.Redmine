#
# Makefile for LArPandoraAlgorithms
#
include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)

LIBCCFILES  := $(wildcard LArClusterAssociation/*.cc) $(wildcard LArClustering/*.cc) $(wildcard LArClusterSeedAssociation/*.cc) $(wildcard LArClusterSplitting/*.cc) $(wildcard LArHelpers/*.cc) $(wildcard LArMonitoring/*.cc) $(wildcard LArObjects/*.cc) $(wildcard LArReclustering/*.cc) $(wildcard LArThreeDSeed/*.cc) $(wildcard LArTwoDSeed/*.cc) $(wildcard LArUtility/*.cc) $(wildcard LArVertex/*.cc)
vpath %.cc LArCalculators LArClusterAssociation LArClustering LArClusterSeedAssociation LArClusterSplitting LArHelpers LArMonitoring LArObjects LArReclustering LArThreeDSeed LArTwoDSeed LArUtility LArVertex

PACKAGE_INCLUDE = $(PACKAGE)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk

override LIBLIBS += $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(PANDORA_LIB) -lPandoraSDK -lPandoraMonitoring

override CPPFLAGS += -DMONITORING=1 -I$(PANDORA_INC) -I$(PANDORA_FQ_DIR)/PandoraSDK/include -I$(PANDORA_FQ_DIR)/PandoraMonitoring/include


