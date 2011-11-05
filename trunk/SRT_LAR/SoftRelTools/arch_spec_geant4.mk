#=======================================================================
#
# arch_spec_geant4.mk:
#
# architecture/site specific makefile fragment for clients of the
# geant4 detector simulation and particle transport libraries.
# This expects to find the variables used by GEANT4 during its build
# which can ge setup with a command like:
#
# % source /usr/local/geant4/4.9.2/env.csh
#
#
# For information about geant4 see http://cern.ch/geant4
#
override CPPFLAGS += -I${G4INCLUDE} 

#
# User interface flags from the environment
#
  override CPPFLAGS += -DG4UI_USE 
  override CPPFLAGS += -DG4UI_USE_TCSH 
  override CPPFLAGS += -DG4UI_USE_XM 
  override CPPFLAGS += -DG4UI_USE_XAW 

  override CPPFLAGS += -DG4VERBOSE 

#
# Visualization options from the environment
#
  override CPPFLAGS += -DG4VIS_USE
  override CPPFLAGS += -DG4VIS_USE_DAWNFILE 
  override CPPFLAGS += -DG4VIS_USE_HEPREPFILE 
  override CPPFLAGS += -DG4VIS_USE_OPENGLX 
  override CPPFLAGS += -DG4VIS_USE_OPENGLXM 
  override CPPFLAGS += -DG4VIS_USE_OPENGL 
  override CPPFLAGS += -DG4VIS_USE_RAYTRACER 
  override CPPFLAGS += -DG4VIS_USE_VRMLFILE 
  override CPPFLAGS += -DG4VIS_USE_ASCIITREE 
  override CPPFLAGS += -DG4VIS_USE_GAGTREE 
  override CPPFLAGS += -DG4VIS_USE_DAWN 
  override CPPFLAGS += -DG4VIS_USE_RAYTRACERX 
  override CPPFLAGS += -DG4VIS_USE_VRML 
  override CPPFLAGS += -DG4VIS_USE_OPENGLX 
  override CPPFLAGS += -DG4VIS_USE_OPENGL 

#
# Probably shouldn't need the CLHEP added to the library list like
# this, but here it is. SRT's arch_spec_CLHEP seems to clash with the
# geant4 configuration of the environment. Go with geant4's
# configuration.
#
G4LIBS = \
-lG3toG4 \
-lG4baryons \
-lG4biasing \
-lG4bosons \
-lG4brep \
-lG4csg \
-lG4cuts \
-lG4decay \
-lG4detector \
-lG4detscorer \
-lG4detutils \
-lG4digits \
-lG4emadjoint \
-lG4emhighenergy \
-lG4emlowenergy \
-lG4empii \
-lG4empolar \
-lG4emstandard \
-lG4emutils \
-lG4error_propagation \
-lG4event \
-lG4FR \
-lG4gdml \
-lG4geombias \
-lG4geomBoolean \
-lG4geomdivision \
-lG4geometrymng \
-lG4geomtext \
-lG4gflash \
-lG4gl2ps \
-lG4globman \
-lG4GMocren \
-lG4graphics_reps \
-lG4had_im_r_matrix \
-lG4had_lll_fis \
-lG4had_mod_man \
-lG4had_mod_util \
-lG4had_muon_nuclear \
-lG4had_neu_hp \
-lG4had_preequ_exciton \
-lG4hadronic_ablation \
-lG4hadronic_abrasion \
-lG4hadronic_bert_cascade \
-lG4hadronic_binary \
-lG4hadronic_body_ci \
-lG4hadronic_coherent_elastic \
-lG4hadronic_crosec_ci \
-lG4hadronic_deex_evaporation \
-lG4hadronic_deex_fermi_breakup \
-lG4hadronic_deex_fission \
-lG4hadronic_deex_gem_evaporation \
-lG4hadronic_deex_handler \
-lG4hadronic_deex_management \
-lG4hadronic_deex_multifragmentation \
-lG4hadronic_deex_photon_evaporation \
-lG4hadronic_deex_util \
-lG4hadronic_em_dissociation \
-lG4hadronic_fragm_ci \
-lG4hadronic_HE \
-lG4hadronic_hetcpp_evaporation \
-lG4hadronic_hetcpp_utils \
-lG4hadronic_incl_cascade \
-lG4hadronic_interface_ci \
-lG4hadronic_iso \
-lG4hadronic_LE \
-lG4hadronic_mgt \
-lG4hadronic_proc_ci \
-lG4hadronic_proc \
-lG4hadronic_qgstring \
-lG4hadronic_qmd \
-lG4hadronic_radioactivedecay \
-lG4hadronic_RPG \
-lG4hadronic_stop \
-lG4hadronic_util \
-lG4hadronic_xsect \
-lG4had_string_diff \
-lG4had_string_frag \
-lG4had_string_man \
-lG4had_theo_max \
-lG4hepnumerics \
-lG4hits \
-lG4intercoms \
-lG4ions \
-lG4leptons \
-lG4magneticfield \
-lG4materials \
-lG4mctruth \
-lG4mesons \
-lG4modeling \
-lG4muons \
-lG4navigation \
-lG4OpenGL \
-lG4optical \
-lG4parameterisation \
-lG4partadj \
-lG4partman \
-lG4partutils \
-lG4phys_builders \
-lG4phys_lists \
-lG4procman \
-lG4RayTracer \
-lG4readout \
-lG4run \
-lG4scoring \
-lG4shortlived \
-lG4specsolids \
-lG4tracking \
-lG4track \
-lG4transportation \
-lG4Tree \
-lG4UIbasic \
-lG4UIcommon \
-lG4UIGAG \
-lG4visHepRep \
-lG4vis_management \
-lG4visXXX \
-lG4volumes \
-lG4VRML \
-lG4xrays \

override LDFLAGS  += -L${G4LIB}/$(G4SYSTEM) -L${CLHEP_DIR}/lib ${G4LIBS}

override LOADLIBES += -L$(G4LIB)/$(G4SYSTEM) $(G4LIBS) -L$(CLHEP_DIR)/lib -lCLHEP
#=======================================================================
