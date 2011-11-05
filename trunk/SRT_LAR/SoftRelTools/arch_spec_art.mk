#=======================================================================
#
#  arch_spec_art.mk:
#
# Define ART-specific build environment
#

#=======================================================================
# Provide a mechanism to install job control files into standard locations
#
ifndef JOB_DIR
  JOB_DIR = $(SRT_PRIVATE_CONTEXT)/job/
endif

ifdef JOBFILES
  JOB_dest = $(foreach v, $(JOBFILES),$(JOB_DIR)$v)	
  SRT_PRODUCTS += $(JOB_dest)
endif

$(filter $(JOB_DIR)%, $(JOB_dest)): $(JOB_DIR)% : %

override CFLAGS += -fPIC
override CXXFLAGS += -fPIC

override CPPFLAGS += -I$(ART_INC)
override CPPFLAGS += -I$(CLHEP_INC)
override CPPFLAGS += -I$(CPPUNIT_INC)
override CPPFLAGS += -I$(LIBSIGCPP_INC)/sigc++-2.0
override CPPFLAGS += -I$(LIBSIGCPP_LIB)/sigc++-2.0/include
override CPPFLAGS += -I$(MESSAGEFACILITY_INC)
override CPPFLAGS += -I$(FHICLCPP_INC)
override CPPFLAGS += -I$(CETLIB_INC)
override CPPFLAGS += -I$(CPP0X_INC)
override CPPFLAGS += -I$(BOOST_INC)

override LDFLAGS  += -L$(ART_LIB)
override LDFLAGS  += -L$(CLHEP_BASE)/lib
override LDFLAGS  += -L$(CPPUNIT_LIB)
override LDFLAGS  += -L$(LIBSIGCPP_LIB)
override LDFLAGS  += -L$(MESSAGEFACILITY_INC)
override LDFLAGS  += -L$(FHICLCPP_INC)
override LDFLAGS  += -L$(CETLIB_INC)
override LDFLAGS  += -L$(BOOST_LIB)

codegen: job 

clean: cleanjob

job: $(JOB_dest) 

#$(foreach v,$(SUBDIRS),$v.py)

$(JOB_dest):
	if [ ! -d $(JOB_DIR) ]; then mkdir $(JOB_DIR); fi
	@echo "<**installing JOB file**> $(@F)"
	$(TRACE)rm -f $@
	$(TRACE)cp $< $@

#%.xml:
#	$(TRACE)$(pass-to-subdirs)

cleanjob:
	if [ ! -d $(JOB_DIR) ]; then mkdir $(JOB_DIR); fi
	$(TRACE)rm -f $@

#=======================================================================
