# ======================================================================
#  larana main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larana/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

project(larpandora)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{ART_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

cet_set_compiler_flags(DIAGS VIGILANT
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
  EXTRA_CXX_FLAGS -std=c++11)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcore v0_01_00 )
find_ups_product( lardata v0_01_00 )
find_ups_product( larevt v0_01_00 )
find_ups_product( larsim v0_01_00 )
find_ups_product( nutools v0_00_00 )
find_ups_product( art v1_06_00 )
find_ups_product( pandora v00_13 )
find_ups_product( cetbuildtools v3_05_00 )

cet_check_gcc()

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

include_directories(LArPandoraAlgorithms)
cet_find_library( PANDORASDK NAMES PandoraSDK PATHS ENV PANDORA_LIB )
cet_find_library( PANDORAMONITORING NAMES PandoraMonitoring PATHS ENV PANDORA_LIB )


# source
add_subdirectory(LArPandoraAlgorithms)
add_subdirectory(LArPandoraInterface)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

