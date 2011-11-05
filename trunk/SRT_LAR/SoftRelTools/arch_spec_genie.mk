#=======================================================================
#
#  arch_spec_genie.mk:
#
# architecture/site specific makefile fragment for clients of GENIE
# neutrino interaction generator interface
#
# For information about GENIE see http://www.genie-mc.org
#
#dont forget the includes and libs for the stuff GENIE depends on
override CPPFLAGS += -I$(LOG4CPP_INC)
override LDFLAGS  += -L$(LOG4CPP_FQ_DIR)/lib -llog4cpp

override CPPFLAGS += -I$(LIBXML2_INC)/libxml -I$(LIBXML2_INC)
override LDFLAGS  += -L$(LIBXML2_FQ_DIR)/lib -lxml2

override CPPFLAGS += -I$(PYTHIA_INC)
override LDFLAGS  += -L$(PYTHIA_FQ_DIR)/lib -lPythia6

override CPPFLAGS += -I$(LHAPDF_INC)/LHAPDF
override LDFLAGS  += -L$(LHAPDF_FQ_DIR)/lib -lLHAPDF

override LOADLIBES += -L$(LOG4CPP_FQ_DIR)/lib -llog4cpp
override LOADLIBES += -L$(LIBXML2_FQ_DIR)/lib -lxml2
override LOADLIBES += -L$(PYTHIA_FQ_DIR)/lib -lPythia6
override LOADLIBES += -L$(LHAPDF_FQ_DIR)/lib -lLHAPDF

GENIELIBS = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

override CPPFLAGS += -I$(GENIE)/src
override LDFLAGS  += $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

override LOADLIBES += $(GENIELIBS)

#$(warning $(LOADLIBES) are libs to be loaded)


#=======================================================================
