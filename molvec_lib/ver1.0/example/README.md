When you run, you first need to mak sure you have the following packages or something equivarent. 

tested environments
|----------------------|
|intel/19.1.2          |
|----------------------|
|openmpi/intel/4.1.1   |
|----------------------|
make

Then before you run make sure you have the right weight file. If you use urea example copy weight_urea.txt as weight.txt. If you use nicotinamide, copy weight_nicotinamide.txt as weight.txt
The following is the sample code for running urea and nicotinamide. 
If the box is orthogonal box... use the keyword lammps_vec
./sa traj_urea.dump lammps_vec INPUT_urea.dat

If the box is non-orthogonal box... use the keyword lammps_vec_tri
./sa nicotinamide_crystal_w_all_info.dump lammps_vec_tri INPUT_nico_crystal.dat

You will get the following output:

atominfo.txt:		providing atoms basic information
new_traj.dump:		if symf_calc_or_NN from INPUT file was NN, it will show NN result for each atom
run_ave_traj.dump:	this result shows running average of NN result. If symf_calc was selected in INPUT file, this file will have only position information.  
symf_result.txt:	if symf_calc_or_NN from INPUT file was symf_calc, this file will store the symmetry function result. 
