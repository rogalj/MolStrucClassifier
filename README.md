# classifierCVs

**[Documentation](https://molstrucclassifier.readthedocs.io/en/latest/)** | **[arXiv](https://arxiv.org/abs/2404.00155)**

<p align="justify"> 
  ClassifierCVs is a local structure classification tool based on NN 
</p>

The following paper describe details of the local structure classification. 

 <link rel="shortcut icon" type="image/png" href="logo.png">


## Table of Contents
- [Introduction](#Introduction)
- [Requirement](#Requirement)
- [Installation](#Installation)
- [Setting up](#SettingUp)
  - [General Set Up](#GeneralSetUp)
  - [Tunable Parameters](#TunableParameters)
- [Data Availability](#DataAvailability)
- [Note](#Note)
- [Author](#Author)
- [Citation](#Citation)

## Introduction
<p align="justify">
  Classifying the local polymorphic structures of molecular materials is a challenging task for both simulation and experiment. We demonstrate a novel approach to for characterizing the local environments in different polymorphs using symmetry function representations. We develop and train a model based on a symmetry function representation, encoding information about the local environment for each molecule. We demonstrate very high classification accuracy for the approach on urea and nicotinamide crystal polymorphs, and extensions to clusters and interfaces. This approach could be applied to new molecules and different complex topologies, providing insights into complex condensed matter phenomena.
  
</p>
 
## Requirement

<p align="justify">
  This package requires C++ compiler. 
  Environments under intel/19.1.2 and openmpi/intel/4.1.1 were tested 
</p>


 

## Installation

<p align="justify">
  To install classifierCVs code, 
  
  ```bash
git clone git@github.com:rogalj/classifierCVs.git
```
Before you run the code make sure you load the c++ environment. 

  ```bash
module load intel/19.1.2 
```

Before you run the code, make sure you compile it by using make.
  ```bash
make 
```

</p>

## Setting Up

### General Set Up
<p align="justify">
  If you already have a dataset with parameters, copy the dataset, neural network weight file, and INPUT file to the working directory. This code accepts a dump file as input trajectory data. In the dump file, you need to have a unique number (sym_type) assigned to atoms, which is used for center of the molecule and point vector representation. For example, urea has 8 atoms. Let's say we use carbon, oxygen, and two nitrogens for the center of the molecule and point vector representations; we can assign 1 for carbon, 2 for oxygen, 3 for the first nitrogen, 4 for the second nitrogen, and 0 for rest of atoms. If you use LAMMPS, we have codes that read/write data file and codes that write dump file with sym_type information. Check the <a class="reference external" href="https://github.com/rogalj/MolStrucClassifier/tree/main/USER-ATOM-VEC-FULL2">USER-ATOM-VEC-FULL2 directory</a> 
  After you build the program using make, you can submit the code as following.

For the orthogonal system, 

```bash
  ./sa traj.dump lammps_sym_new INPUT.dat
```

For the non-orthogonal system, 

```bash
  ./sa traj.dump lammps_sym_new_tri INPUT.dat
```


If you use <a class="reference external" href="https://github.com/rogalj/classifierCVs/tree/main/molvec_lib/ver2.0">ver2.0</a>, you can also run the code without sym_type if the atom type you chose for the center molecule and point-vector representation are unique numbers. In this case You can run the following code. 

For the orthogonal system, 

```bash
  ./sa traj.dump lammps_vec INPUT.dat
```

For the non-orthogonal system, 

```bash
  ./sa traj.dump lammps_vec_tri INPUT.dat
```
</p>

### Tunable parameters
The following is the first part of the tunable parameter that stores the general information. The following is the example using nicotinamide. 
<p align="justify">

```
num_sf           24                     //This should be the sum of all symmetry function
num_pvsf_type    2                      //number of point vector symmetry function
triclinic        true                   //true if tye system is triclinic. false if it is not.
center atom      2                      //sym_type of the center atom. If you are not using sym_type, this should be the atom type This number need to be unique
atm in molec     16                     //number of atoms in the molecule If you have more than one morecule type, choose the higher number
cutoffmin        6.8                    //cutoff rmin
cutoffmax        7.0                    //cutoff rmax
symf_calc_or_NN  NN                     //symf_calc: calculate only symmetry function and print out in symf_result.txt. This can be used for the neural network input
                                        //NN: calculate NN using the existing weight file.
NNrdf            true                   //Calculate the distance from surface. true for yes, false for no
fullNeigh        7                      //Number of neighbors the form usually get. You can adjust this number to change the surface info.  
rdiff            0.2                    //default value is 0.2
```

The following is the second part of the tunable parameter that stores detailed symmetry function information. Consider each 5-6 lines as block. number_of_para represents the number of the symmetry function of the same type. "point" is a parameter that indicates whether the function uses point representation or point vector representation. If you choose point representation, the subsequent section should be "G2G3". If you choose pointvector representation, the subsequent section should be "vec_id" which stores sym_type or atom type you would like to use for vector. In "G2G3", choose which symmetry function you would like to use. For "Rs", "eta", and, "kappa" valus write the same number of parameter as you defined in number_of_para. 

```
number_of_para   4                      
point            true
G2G3             G2
Rs               3.75 4.86 4.9 5.9
eta              1.26 2.11 0.1 0.016

number_of_para   4
point            false
vec_id           2 1
G2G3             G2
Rs               0.01 2.71 2.9 1.69
eta              0.66 1.11 0.9 4.0

number_of_para   4
point            false
vec_id           3 4
G2G3             G2
Rs               0.56 2.26 1.18 1.15
eta              3.0 1.11 17.2 2.2

number_of_para   4
point            true
G2G3             G3
kappa            1.06 0.51 1.03 2.41

number_of_para   4
point            false
vec_id           2 1
G2G3             G3
kappa            0.13 7.91 5.05 3.2


number_of_para   4
point            false
vec_id           3 4
G2G3             G3
kappa            1.05 2.21 2.36 13.22
```

</p>

## Note

<p align="justify">
We use NYU high performance computer GREENE as an example. If you are NYU student/employee who does not have access to the account, go to `NYU HPC account page <a class="reference external" href="https://sites.google.com/nyu.edu/nyu-hpc/accessing-hpc/getting-and-renewing-an-account?authuser=0">NYU HPC account page</a>  to get access. 
</p>


## Author
 
* Jutta Rogal
  * New York University / Freie Universität
  * jutta.rogal@nyu.edu

* Daisuke Kuroshima
  * New York University 
  * daisuke.Kuroshima@nyu.edu
 

