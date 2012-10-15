#
# $Id: GNUmakefile,v 1.4 2009/06/05 18:49:52 t962cvs Exp $
#
# messier@indiana.edu
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := BeamAna
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE) 


########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lUtilities -lRawData -lT962_Objects
