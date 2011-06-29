#
# $Id: GNUmakefile,v 1.3 2008/08/26 21:20:33 brebel Exp $
#
# brebel@fnal.gov
#
include SoftRelTools/arch_spec_root.mk

SUBDIRS     := T962_Objects MergeData VertexActivity MatchTracks MergeScan ScanFilter Kinematics SimKinematics T962G4Ana AnalysisTree TGMuon CCQEanalysis MuonAna

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service
