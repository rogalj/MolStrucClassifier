# Example files

<p align="justify"> 


## input file
Here, we have nicotinamide as sample. You need to have the run file, data file, init file, settings file to execut the program. Make sure you compiled atom.cpp, atom.h, atom_vec_full2.cpp, atom_vec_full2.h, dump_custom.cpp, dump_custom.h, dump_custom_mpiio.cpp, and dump_custom_mpiio.h in the lammps' src directory. 
In the data file you need to have the following format in the Atoms section in the data file. 

```
Atoms  # full2

atom ID, molID, atom_type, charge, sym_type, x, y, z
```
This is the sample line for Atoms information. 
```
1 1 5 -0.326615 1 2.785 7.001 1.945
```
In the system.init file, you need to change the atom_style to full2. (first line in the file)

**run_G2G3nico.lmp**: This file contains the simulation detail. 

**system.data**: This data contains atoms information such as coordinate, box size, masses, charges etc. 

**system.init**: This file contains the function detail of the system

**system.settings**L This file contains force field information.

**DK_job.slurm**: slurm file that calls LAMMPS. 


> [!NOTE]
> When you make the sym_type, you need to have a unique number for the one you use for the symmetry function calculation. For others, you can use the same number. In this example, I used 0 for atoms that I did not use for the calculation.


## output file
After the calculation, you will get the following files. 

**traj.dump**: This file contains trajectory information of the simulation. Each sym_type is saved before coordinate information. 

**screen.log** and **log.lammps**: These files showed the summary of the LAMMPS simulation. 
