# ver 2.0 is now available 

<p align="justify"> 

  Here we have codes that calculate the symmetry function and NN by itself or using LAMMPS. 
  

## Preparing data

### Trajectory file
If you wish to get the symmetry function result or neural network result using symmetry function, you can use this code. Currently this code accepts a dump file as input. 
Inside the dump file, you need to have the data in the following order. 
  
  ITEM: ATOMS id mol element x y z vx vy vz

where the element needs to be a number and elements that will be used for the symmetry function need to be unique in the molecule. (If you use <a class="reference external" href="https://github.com/rogalj/classifierCVs/tree/main/molvec_lib/ver2.0">ver2.0</a>, you do not need to worry about this)

### input file 

input file stores all the tunable parameters. The detail of the inputs are shown <a class="reference external" href="https://github.com/rogalj/classifierCVs/tree/main?tab=readme-ov-file#TunableParameters">here</a>.

### weight file 

The weight file is necessary in the running directory. If you wish to calculate the symmetry function without neural network, you can use the sample weight file instead. Make sure you change the first number and last number of the second line to number of input and number of output. 


## Stand Alone


Firstly, you need to load the module as following

  ```bash
module load openmpi/intel/4.1.1
```

Then you need to load files using 

```bash
make
```


When you run the code you need to have the following command. 

  ```bash
./sa trajectory_file file_type  name_of_input_data.data
```

The file_type need to be one of the following. 


| file_type     | Description |
| ------------- | ------------- |
| lammps_vec  | If you have orthogonal box, use this command |
| lammps_vec_tri  | If you have triclinic box box, use this command |












## LAMMPS
  
  If you would like to use this with LAMMPS, make sure you have the following LAMMPS extentions.  
  
  
   |Package        |
   |---------------|
   | KSPACE        |
   | MANYBODY      |
   | MOLECULE      |
   | OPT           |
   |EXTRA-MOLECULE |

  To connect this libraly to LAMMPS, you need to make the soft links. Inside of the libraly(this) directory, make DAFED.inc file that contains the following lines.
```bash
  DAFED_LOAD= /your_libraly_directory/libcv_nn.so -ldl
  DAFED_DEPENDENCIES = /your_libraly_directory/llibcv_nn.so
```

You can then load the codes using the following commands.
```bash
make lib
```

Then, inside the lammps directory (above src direcotry), type the following commands. 

```bash
  ln -s /your_libraly_directory/Dafed.inc Dafed.inc
  ln -s /your_libraly_directory/dafed_link.h Dafed.h
```

Nextly, you need to go to src directory in LAMMPS and open Makefile.package. Then add the following line
```bash
PKG_LIB =    $(DAFED_LOAD)
```
Lastly, you need to go to MAKE directory in the src directory. and replace the line 10 to 
```bash
CCFLAGS =       -g -O3 -std=c++11 -restrict
```  
  





</p>
