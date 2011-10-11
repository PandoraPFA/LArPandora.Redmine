export GROUP=argoneut
. /grid/fermiapp/common/tools/setup_condor.sh
# to run one copy on grid.
#jobsub -g -q condor_uBdetMC.sh `whoami`  `pwd`
# -N 197/198 for all ArgoNeuT runs
# -N 33 for neutrino-mode
jobsub -X509_USER_PROXY  /scratch/kpartyka/grid/kpartyka.argoneut.proxy -g -N 33 -dOUT /argoneut/data/users/kpartyka/out -q condor_datakinreco_grid.sh `whoami`  `pwd`

# Then, condor_q USER to see the status.
