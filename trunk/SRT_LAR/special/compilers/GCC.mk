# Use the g++ command to build shared libraries
# (to include the additional C++-specific information)
override SHAREDAR := g++

# Tips:
# CXXFLAGS - used when compiling C++ programs
# CCFLAGS - used when compiling C programs
# LDFLAGS - used when linking binary programs
# ARFLAGS - used when linking static libraries
# SHAREDARFLAGS - used when linking shared libraries

# See $SRT_PUBLIC_CONTEXT/SoftRelTools/include/arch_spec.mk
# for more flags/options.

# For LArSoft, CXXFLAGS and SHAREDARFLAGS are the
# most relevant.

# Compile using the debug flag. Turn off
# optimization, which can sometimes interfere
# with line-by-line debugging.
ifeq ($(findstring debug,$(SRT_QUAL)),debug)
   override CXXFLAGS += -g -O0 -Wall -Wwrite-strings -Wno-inline -Woverloaded-virtual
   override CCFLAGS  += -g -O0 -Wall -Wwrite-strings -Wno-inline -Woverloaded-virtual
   override LDFLAGS  += -g
   override ARFLAGS  += -g
   override SHAREDARFLAGS  += -g
endif

# If SRT_QUAL contains "prof", compile with 
# both the debugger and the profiler options.
ifeq ($(findstring prof,$(SRT_QUAL)),prof)
  override CXXFLAGS += -g -pg -O0 -Wall -Wwrite-strings -Wno-inline -Woverloaded-virtual
  override CCFLAGS  += -g -pg -O0 -Wall -Wwrite-strings -Wno-inline -Woverloaded-virtual
  override LDFLAGS  += -g -pg 
  override ARFLAGS  += -g -pg
  override SHAREDARFLAGS  += -g -pg
endif
