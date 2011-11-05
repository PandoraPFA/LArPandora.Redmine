#=======================================================================
#
#  arch_spec_xercesc.mk:
#
#  architecture/site specific makefile fragment for clients of
#  Xerces-c XML interface
#

override CPPFLAGS += -I${XERCESCROOT}/include
override LDFLAGS  += -L${XERCESCROOT}/lib -lxerces-c

override LOADLIBES += -L${XERCESCROOT}/lib -lxerces-c

#=======================================================================
