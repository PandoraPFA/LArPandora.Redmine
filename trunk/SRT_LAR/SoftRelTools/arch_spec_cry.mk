#=======================================================================
#
#  arch_spec_cry.mk:
#
# architecture/site specific makefile fragment for clients of CRY
# cosmic ray generator interface
#
# For information about CRY see http://nuclear.llnl.gov/simulation
#

override CPPFLAGS += -I${CRYHOME}/src
override LDFLAGS  += -L${CRYHOME}/lib -lCRY

override LOADLIBES += -L$(CRYHOME)/lib -lCRY

#=======================================================================
