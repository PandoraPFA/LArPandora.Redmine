###
#
#  arch_spec_root.mk:
#
#      architecture/site specific makefile fragment for clients of ROOT.
#
# This WILL only build in the shared library directory (see cintworkdir)

ROOTLIBS  = $(shell root-config --libs) -lEG -lTreePlayer -lGeom -lFFTW -lReflex -lEGPythia6

ROOTGLIBS =  -lGui

override CPPFLAGS += -DUSE_ROOT

ifeq (,$(ROOTSYS))  
  override CPPFLAGS += $(shell root-config --cflags) 
  override LDFLAGS += $(shell root-config --libs) 
else 
  override CPPFLAGS += -I$(ROOTSYS)/include 
  override LDFLAGS += -L$(ROOTSYS)/lib 
endif

override LOADLIBES += $(ROOTLIBS)
ifndef USE_ROOT_NO_GRAPH
override LOADLIBES += $(ROOTGLIBS)
endif

ifneq (,$(findstring IRIX6,$(BFARCH)))
    override LOADLIBES += -lgen
    override ROOTLIBS  += -lgen
endif

#
#  CINT dictionary (create it with *Cint.cc extention not to mix with the rest 
#  of the sources) 
#
ifdef ROOTCINT

ifeq (static,$(LIB_TYPE))
cintworkdir = $(staticlib_o_dir)
else
cintworkdir = $(sharedlib_o_dir)
endif

cintdir = $(findstring test,$(shell pwd | sed -e 's%$(SRT_PRIVATE_CONTEXT)/%%'))

CINTOUTPUT = $(PACKAGE)$(cintdir)Cint.cc
CINTOBJECT = $(PACKAGE)$(cintdir)Cint.o
CINT_DEPENDS       = $(cintworkdir)$(PACKAGE)$(cintdir)Cint.d
CINT_DEPENDS_1     = $(cintworkdir)$(PACKAGE)$(cintdir)Cint.ccd

# Generate automatically the list of header that should be parsed
# by cint.
CINT_CXXFILES = $(wildcard *.cxx)
CINTLIST += $(addsuffix .h, $(basename $(CINT_CXXFILES)))

# Tell make that $(CINTOUTPUT) can be found in cintworkdir
vpath $(CINTOUTPUT) $(cintworkdir)

# Now add the object to the list of object files for the 
# library.  If SoftRelTools internals change this may need to 
# change
actual_sharedlib_files += $(cintworkdir)$(CINTOBJECT)
$(SHAREDLIB): $(cintworkdir)$(CINTOBJECT)

actual_staticlib_files += $(cintworkdir)$(CINTOBJECT)
$(STATICLIB): $(cintworkdir)$(CINTOBJECT)

ifndef NO_GENERATE_DEPENDS
depend: $(CINT_DEPENDS) $(CINT_DEPENDS_1)
endif #NO_GENERATE_DEPENDS

codegen: $(cintworkdir)$(CINTOUTPUT)

$(cintworkdir)$(CINTOUTPUT): $(CINTLIST) LinkDef.h
	echo "Generating $(@F) dictionary..."
	$(check_dep_dir)
	rootcint -f $@ -c -p \
		$(CPPFLAGS) $(CINTLIST) LinkDef.h


define ccd_generate_depends
	grep -v '.*Cint.cc' $(dir $@)/$(basename $(notdir $<)).d | sed -e 's?Cint.o?Cint.cc?g' > $(dir $@)/$(basename $(notdir $<)).ccd \
	|| /bin/rm -f $(dir $@)/$(basename $(notdir $<)).ccd
endef

#
# The following defines rules to create the dependency files
# for source files that are generated (by codegen for examples).
# this is no root specific and propably belongs into a
# post_standard.mk file.  The only restriction (now) is
# that the file have a suffix of Cint.cxx
#

libccdepends := $(patsubst %.cc,%.ccd, $(wildcard $(cintworkdir)$(PACKAGE)$(cintdir)Cint.cc))

ifneq ("$(libccdepends)","")
    -include $(libccdepends)
endif

ifndef NO_GENERATE_DEPENDS

# We have this rule because we NEED this .d file to create the .ccd
# EVEN when the file has not been compiled yet. 
$(cintworkdir)$(PACKAGE)$(cintdir)Cint.d: $(cintworkdir)$(PACKAGE)$(cintdir)Cint.cc
	echo "<**special depend**> $(<F)"
	$(check_dep_dir)
	$(cxx_generate_depends)

$(cintworkdir)$(PACKAGE)$(cintdir)Cint.o: $(cintworkdir)$(PACKAGE)$(cintdir)Cint.cc

$(cintworkdir)$(PACKAGE)$(cintdir)Cint.ccd: $(cintworkdir)$(PACKAGE)$(cintdir)Cint.d
	echo "<**depend**> $(@F)"
	$(check_dep_dir)
	$(ccd_generate_depends)

endif #NO_GENERATE_DEPENDS

endif # ROOTCINT
