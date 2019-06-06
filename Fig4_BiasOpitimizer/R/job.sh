#!/bin/bash -l
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -l walltime=12:00:00
#PBS -o log/SWIFT.o$PBS_JOBID
#PBS -e log/SWIFT.e$PBS_JOBID

ml R/3.4.4-intel-2018a-X11-20180131; ulimit -s unlimited
cd /kyukon/data/gent/vo/000/gvo00074/felicien/R/InvSWIFToptim/Fig4_BiasOpitimizer

echo "source('./R/analsysisOptimization.R')" | R --vanilla

