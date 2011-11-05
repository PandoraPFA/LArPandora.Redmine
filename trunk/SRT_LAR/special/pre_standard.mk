#=======================================================================
#
#  arch_spec_art.mk:
#
# Define ART-specific build environment
#

#=======================================================================
#if we have the genreflex source file, we add what it makes to 
# the "codegen" target
#=======================================================================

ifneq (,$(wildcard classes_def.xml))

  MWM_CLASSES_DEF_GEN= $(PACKAGE)_dict.cc  $(PACKAGE)_map.cc
  codegen: $(MWM_CLASSES_DEF_GEN) $(warning genreflex case found)

   $(MWM_CLASSES_DEF_GEN): classes_def.xml
	genreflex classes.h -s classes_def.xml -o $(PACKAGE)_dict.cc $(CPPFLAGS) --deep --fail_on_warnings --iocomments  --capabilities=classes_ids.cc  -D_REENTRANT -DGNU_SOURCE  -DGNU_GCC  && mv classes_ids.cc $(PACKAGE)_map.cc

#   lib:$(shlibdir)lib$(PACKAGE)_dict_plugin.so $(shlibdir)lib$(PACKAGE)_map_plugin.so

#   $(shlibdir)lib$(PACKAGE)_dict_plugin.so: $(PACKAGE)_dict_plugin.o
#	echo "<**linking dict_plugin**>"
#	cd $(sharedlib_o_dir); $(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $(shlibdir)lib$(PACKAGE)_dict_plugin.so $(PACKAGE)_dict_plugin.o $(LIBLIBS)
#   
#   $(shlibdir)lib$(PACKAGE)_map_plugin.so: $(PACKAGE)_map_plugin.o
#	echo "<**linking map_plugin**>"
#	cd $(sharedlib_o_dir); $(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $(shlibdir)lib$(PACKAGE)_map_plugin.so $(PACKAGE)_map_plugin.o $(LIBLIBS)

    SRT_PRODUCTS += $(MWM_CLASSES_DEF_GEN)
endif

#=======================================================================
#if we have _plugin  source files, we add an extra .so to the dependencies of the lib target.
#=======================================================================

PLUGINSRC=$(sort $(wildcard *_dict.cc *_map.cc *_module.cc *_service.cc *_source.cc) $(MWM_CLASSES_DEF_GEN))

ifneq (,$(PLUGINSRC))
  PLUGINOBJ=$(addprefix $(shlibdir),$(addsuffix .o,$(basename $(PLUGINSRC))))
  PLUGINLIB=$(addprefix $(shlibdir)lib,$(addsuffix .so,$(basename $(PLUGINSRC))))
  SRT_PRODUCTS += $(PLUGINLIB) #$(warning Plugin library base $(PLUGINSRC) for $(PLUGINLIB) found)
#  $(warning $(LIBLIBS) is LIBLIBS)

  lib%_dict.so : %_dict.o
	echo "<**linking**> $(@F)"
	$(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $@ $<  $(LIBLINK) $(LIBLIBS)

  lib%_map.so : %_map.o
	echo "<**linking**> $(@F)"
	$(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $@ $<  $(LIBLINK) $(LIBLIBS)

  lib%_module.so : %_module.o
	echo "<**linking**> $(@F)"
	$(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $@ $<  $(LIBLINK) $(LIBLIBS)

  lib%_service.so : %_service.o
	echo "<**linking**> $(@F)"
	$(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $@ $<  $(LIBLINK) $(LIBLIBS)

  lib%_source.so : %_source.o
	echo "<**linking**> $(@F)"
	$(SHAREDAR) $(SHAREDARFLAGS) $(SHAREDAROFLAG) $@ $<  $(LIBLINK) $(LIBLIBS)

  lib: $(PLUGINLIB) #$(warning Saw plugin sources.. $(PLUGINSRC))

  all: lib

  lib: $(sharedlib_o_dir)
	mkdir -p $(sharedlib_o_dir)

endif


#=======================================================================
