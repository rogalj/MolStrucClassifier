#!/bin/bash

#SBATCH --job-name=t31-100-pbe
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
module purge
module load intel/19.1.2
module load openmpi/intel/4.1.1

echo "Start"
date
hostname
cd $SLURM_SUBMIT_DIR
pwd
#lmp_test_Dec -in run_symftype.lmp -screen screen.log
#lammps -in run_symftype_w_boundary.lmp -screen screen.log
#lammps -in run_symftype_w_o_boundary.lmp -screen screen.log
#lammps -in run_NNout_w_boundary.lmp -screen screen.log
#lammps -in run_NNout_w_o_boundary.lmp -screen screen.log
#lammps -in run_G2G3.lmp -screen screen.log
#lammps -in run_NNout.lmp -screen screen.log
#lammps -in run_G2G3.lmp_original -screen screen.log
#lammps -in run_NNout.lmp -screen screen.log
#lammps -in run_symftype.lmp -screen screen.log
#mpirun -np 8 lammps -in run_NNout.lmp -screen screen.log
srun lammps -in run_NNout.lmp -screen screen.log
#srun lammps -restart2data restart_from0.13000000.% restart_file_for_weight.data

echo "Finished"
date
### run your executable progra
######################
