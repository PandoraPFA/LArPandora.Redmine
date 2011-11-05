#!/bin/bash
#
# Author: echurch@fnal.gov from dbox@fnal.gov
# .fcl-specific ART codework by joshua.spitz@yale.edu
# A script to run argoneut framework (ART) jobs on the local cluster/grid.
#
# It takes 4 arguments, the latter two of which are specified. First two are inherited.:
#   cluster    - the condor job cluster number
#   process    - the condor job process number, within the cluster
#   user       - the username of the person who submitted the job
#   submitdir  - the directory from which the job was submitted ( not used at this time).
#
# Outputs:
#  - All output files are created in the grid scratch space.  At the end of the job
#    all files in this directory will be copied to:
#      /grid/data/argoneut/outstage/$user/${cluster}_${process}
#    This includes a copy of the input files.
#
# Notes:
#
# 1) For documentation on using the grid, see
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtm
#    For details on cpn and outstage see:
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml#cpn
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml#outstage
#
sleep `echo $((RANDOM%1000+0))`

verbose=T

# Copy arguments into meaningful names.
cluster=${CLUSTER}
process=${PROCESS}
user=$1
submitdir=$2
echo "Input arguments:"
echo "Cluster:    " $cluster
echo "Process:    " $process
echo "User:       " $user
echo "SubmitDir:  " $submitdir
echo " "

# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
# Next 3 commented-out lines are outdated cruft for grid, now replaced at 
# Dennis's suggestion. EC, 23-Nov-2010. 
ORIGDIR=`pwd`
#TMP=`mktemp -d ${OSG_WN_TMP:-/var/tmp}/working_dir.XXXXXXXXXX`
#TMP=${TMP:-${OSG_WN_TMP:-/var/tmp}/working_dir.$$}
TMP=`mktemp -d ${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0

# End of the section you should not change.

# Directory in which to put the output files.
outstage=/argoneut/data/outstage/$user


cd /argoneut/app/users/saima/lar_sim/sam_batch

touch /argoneut/data/outstage/saima/simkinlist
chmod g+rw /argoneut/data/outstage/saima/simkinlist

# if [ "$i" ]
# then
# i=$(( $i + 100 ))
# else
# i=0
# fi


# ls -lrt /argoneut/data/outstage/spitz7/469154*/mergedata.root|awk '{print $5,$9}' >> ~/list
# sort ~/list|uniq -w32 -d |awk '{print $2}' >> ~/dupelist
# cat ~/dupelist | while read line; do rm -f "${line}"; done
# rm ~/list
# rm ~/dupelist



for i in `ls -d /argoneut/data/outstage/saima/878844_*/mergesim.root|sort`; 
do 
inlist=`grep $i /argoneut/data/outstage/saima/simkinlist`

if [ "$inlist" == "" ]
then

echo $i>>/argoneut/data/outstage/saima/simkinlist




break

fi


done



# rm scanmerge.fcl.dummy

cd $TMP
cp -r $ORIGDIR/* .



# This is needed because of pollution from setup vdt. -- Dennis Box's comment.
# I, however, don't wanna blow out part of LIB that tells ART
# where to pick up modules. EC, 22-Nov-2010.
unset LD_LIBRARY_PATH

# Establish environment and run the job.
#source /grid/fermiapp/products/mu2e/setupmu2e.sh
#source /grid/fermiapp/mu2e/Offline/v0_2_6/setup.sh
export GROUP=argoneut
export EXPERIMENT=argoneut
export HOME=/argoneut/app/users/saima/lar_sim/sam_batch
## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.
source /grid/fermiapp/lbne/lar/code/larsoft/releases/development/setup/setup_larsoft_fnal.sh
cd /argoneut/app/users/saima/lar_sim; srt_setup -a
cd $TMP
chmod -R g+rw *
# fw thisjob.py >& thisjob.log

# lar -c job/simkinematics_reco.fcl -s /argoneut/app/users/spitz7/larsoft_ART9/t962gana.root --nskip "$i" -n 100  >&thisjob.log
lar -c job/kinematics_reco.fcl -s "$i"  >& thisjob4.log

# Make sure the user's output staging area exists.
test -e $outstage || mkdir $outstage
if [ ! -d $outstage ];then
   echo "File exists but is not a directory."
   outstage=/grid/data/argoneut/outstage/nobody
   echo "Changing outstage directory to: " $outstage 
   exit
fi
chmod -R g+rw *
# Make a directory in the outstage area to hold all files from this job.
mkdir ${outstage}/${cluster}_${process}
chmod -R g+rw *
stat=`echo $?`

if [ "$stat" == "0" ]
then
/grid/fermiapp/minos/scripts/cpn * ${outstage}/${cluster}_${process}
fi
chmod -R g+rw $outstage/${cluster}_${process}


exit 0
