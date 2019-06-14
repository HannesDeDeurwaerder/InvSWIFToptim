#!/bin/bash -l
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -l walltime=1:00:00
#PBS -o R/log/SWIFT.o$PBS_JOBID
#PBS -e R/log/SWIFT.e$PBS_JOBID

ml R; ulimit -s unlimited
cd /kyukon/home/user/gent/408/vsc40883/InvSWIFToptim/Fig4_BiasOpitimizer

echo "source('./R/run_parallel.r')" | R --save

