# arch_spec_CLHEP.mk
#
# CLHEP_DIR, CLHEP_LIB and CLHEP_INC are environment variables
#   containing (or optionally loaded with) the
#   location of the CLHEP files on the local machine.
#
# Original version:
# Lars Rasmussen, 16-Sep-98
#

extpkg:=CLHEP
CLHEP_DIR_DEFAULT := /usr/local
ifndef CLHEP_BASE_DIR
    arch_spec_warning:=\
    "Using default value CLHEP_DIR = $(CLHEP_DIR_DEFAULT)"
    CLHEP_BASE_DIR := $(CLHEP_DIR_DEFAULT)
endif
ifndef CLHEP_INC 
    CLHEP_INC = $(CLHEP_BASE_DIR)/include
endif
ifndef CLHEP_LIB
    CLHEP_LIB = $(CLHEP_BASE_DIR)/lib
endif

include SoftRelTools/specialize_arch_spec.mk

#override CPPFLAGS  += -I$(CLHEP_INC)
#override LDFLAGS   += -L$(CLHEP_BASE)/lib
override LOADLIBES += -lCLHEP
