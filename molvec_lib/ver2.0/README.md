# Ver2.0
<p align="justify"> 


  
In this version, we use sym_type as a input for the symetry function. This makes more flexibility on input molecules. You can still use the ver1.0 code in here as well. Make sure you use the right file type

If you use this, make sure you have sym_type in trajectory file as following.

OLD

ITEM: ATOMS id mol element x y z vx vy vz

NEW

ITEM: ATOMS id mol element ${\textsf{\color{red} sym\\_type}}$  x y z vx vy vz

> [!NOTE]
> If you want to create a dump file using LAMMPS, you can use the following LAMMPS extension to insert sym_type
> 
> <a class="reference external" href="https://github.com/rogalj/classifierCVs/tree/main/USER-ATOM-VEC-FULL2">USER-ATOM-VEC-FULL2</a> into the data file.


After you load the module, as following
sample command
  ```bash
module load openmpi/intel/4.1.1
```
You have two options for the new style. Keep in mind that you can still use the old version too. 

  ```bash
./sa_test trajectory_file file_type  name_of_input_data.data
```

List of available file type

| file_type     | Description |
| ------------- | ------------- |
| lammps_vec  | If you use old version with orthogonal box, use this command |
| lammps_vec_tri  | If you use old version with triclinic box box, use this command |
|lammps_sym_new | NEW If you use new version with orthogonal box, use this command|
|lammps_sym_new_tri|NEW If you use new version with triclinic box, use this command|


</p>
