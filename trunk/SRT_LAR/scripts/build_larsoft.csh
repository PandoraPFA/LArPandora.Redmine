#!/bin/csh

# Wrapper Script to build the LArSoft software.
# This sources the development release of the software
# and updates all of the packages.
# It then calls a perl script to make a clean build 
# in standard and debug modes.
# Stolen from NOvA, script originally written
# February 12 2011 - Gavin S. Davies
#
# brebel@fnal.gov

#set the HOSTNAME explicitly because it doesn't seem to get picked up otherwise
#it is ok to set it like this because we will always build the fnal version on lbnegpvm01.fnal.gov
setenv HOSTNAME lbnegpvm01.fnal.gov

source /grid/fermiapp/lbne/lar/code/larsoft/setup/setup_larsoft_fnal.csh

echo update $SRT_DIST/setup
cd $SRT_DIST/setup
svn update

set logfile="/grid/fermiapp/lbne/lar/code/larsoft_make.log"
set debuglogfile="/grid/fermiapp/lbne/lar/code/larsoft_debug_make.log"

echo "Starting FNAL LArSoft Update at " > $logfile
date >> $logfile 
$SRT_PUBLIC_CONTEXT/SRT_LAR/scripts/lar_update_rel -rel development >> $logfile 
echo "Output written to" $logfile
cat $logfile

#do a verbose build to catch all warnings
setenv VERBOSE 1

$SRT_PUBLIC_CONTEXT/SRT_LAR/scripts/lar_build >> $logfile
echo done at `date`

echo "Building in debug mode"
echo "Debug output written to" $debuglogfile

#the following 2 lines are ugly hacks because when run as a cronjob this script
#croaks on srt_setup SRT_QUAL=debug for some reason
setenv SRT_QUAL debug
setenv SRT_SUBDIR Linux2.6-GCC-debug

$SRT_PUBLIC_CONTEXT/SRT_LAR/scripts/lar_build -debug >& $debuglogfile
echo done at `date`

$SRT_PUBLIC_CONTEXT/SRT_LAR/scripts/make_lar_build_log $debuglogfile

