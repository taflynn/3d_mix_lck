#!/bin/bash
#SBATCH -c 8
#
# number of nodes
#SBATCH -N 1
#
# node
#SBATCH -p defq
#
# set the $OMP_NUM_THREADS variable
ompthreads=$SLURM_JOB_CPUS_PER_NODE
export OMP_NUM_THREADS=$ompthreads
#

./gp_lck
