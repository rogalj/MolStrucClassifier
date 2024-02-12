#!/bin/bash

#SBATCH --job-name=nn_fit

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

### Load required modules
module purge
module load openmpi/intel/4.1.1



### set directory for job execution
cd $SLURM_SUBMIT_DIR
pwd


### run programme
echo 'Starting NN fit'
date
time ./ml-classification.exe input.dat > screen.out  
date
echo 'Finished'

######################


