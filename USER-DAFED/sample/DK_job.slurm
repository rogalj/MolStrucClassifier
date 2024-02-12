#!/bin/bash

#SBATCH --job-name=t31-100-pbe
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=0:20:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
##SBATCH --gres=gpu:mi50:1

module purge
module load intel/19.1.2
module load openmpi/intel/4.1.1

echo "Start"
date
hostname
cd $SLURM_SUBMIT_DIR
pwd
srun lammps -in run_DAFED.lmp -screen screen.log

echo "Finished"
date
### run your executable progra
######################
