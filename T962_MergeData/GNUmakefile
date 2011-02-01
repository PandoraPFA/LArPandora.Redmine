#
# $Id: GNUmakefile,v 1.1.1.1 2009/12/07 03:50:52 kpartyka Exp $
#
# brebel@fnal.gov
#
include SoftRelTools/arch_spec_root.mk

LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK    := $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lRawData -l$(PACKAGE)
########################################################################
include SoftRelTools/standard.mk

include SoftRelTools/arch_spec_art.mk

override LIBLIBS += -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service





