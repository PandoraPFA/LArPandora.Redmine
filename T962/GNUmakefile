#
# $Id: GNUmakefile,v 1.3 2008/08/26 21:20:33 brebel Exp $
#
# brebel@fnal.gov
#
include SoftRelTools/arch_spec_root.mk

SUBDIRS     := VertexActivity

LIBLINK    := $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
