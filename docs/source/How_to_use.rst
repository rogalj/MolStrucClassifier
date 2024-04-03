How to use
==========

Symmetry function codes
-----------------------


ver1.0
^^^^^^

Here, we have codes that calculate the symmetry function and NN by itself or using LAMMPS. `ver2.0 <ver2.0_>`_ is available, so make sure you use the right one. 

Data Preparation
""""""""""""""""

To run this code, you need the following files: 

1. trajectory file (.dump)

If you wish to get the symmetry function result or neural network result using the symmetry function, you can use this code. Currently, this code accepts a dump file as input. Inside the dump file, you need to have the data in the following order, which is the same as LAMMPS full atom style. 

.. code-block:: console

   ITEM: ATOMS id mol element x y z vx vy vz

.. warning:: 

   Make sure that the element type is written in numbers not in character or string! Also make sure that the element type you use for symmetry function is unique. If it is not, you cannot use ver1.0. Check `ver2.0 <ver2.0_>`_ instead. 

   
2. Input file

The INPUT.dat file stores all the parameters such as the number of symmetry functions, parameters for the symmetry function, system information, etc. Make sure you check each line and adjust values accordingly. An example of this file is stored in `example <https://github.com/rogalj/MolStrucClassifier/tree/main/molvec_lib/ver1.0/example>`_. 

3. Weight file (only for NN calculation)

The weight file is necessary in the running directory. If you wish to calculate the symmetry function without using the neural network, you can use the sample weight file instead. Make sure you change the first number and last number of the second line to the number of input and the number of output.

Stand Alone
"""""""""""

Firstly, you need to load the module as follows before running the code:

.. code-block:: console

   module load intel/19.1.2

Then, you need to compile files using:

.. code-block:: console

   make

When you run the code you need to have the following command:

.. code-block:: console

   ./sa trajectory_file file_type  name_of_input_data.data

The file_type needs to be one of the following:

.. list-table:: 
   :widths: 25 70
   :header-rows: 1

   * - File Type
     - Description
   * - lammps_vec
     - Use this command if you have an orthogonal box.
   * - lammps_vec_tri
     - Use this command if you have a non-orthogonal box.


Connecting with LAMMPS
""""""""""""""""""""""

If you want to connect this to LAMMPS, please refer to the instructions on how to install and use it in the `Install LAMMPS code <file:///home/drk354/Downloads/MolStrucClassifier-main/docs/build/html/Install.html#install-lammps-code>`_.




Parameters
""""""""""

The following is the basic information you need for the symmetry function code. 

.. list-table::
   :widths: 25 70
   :header-rows: 1

   * - Name 
     - Description
   * - num_sf
     - Number of total symmetry functions
   * - num_pvsf_type
     - Number of point-vector symmetry functions
   * - triclinic
     - True if this is in a non-orthogonal box
   * - center atom
     - Element id of the center atom. Must be unique in the molecule
   * - atm in molec
     - Number of atoms in the molecule. If you have defect, choose the one with more atoms 
   * - cutoffmin
     - Minimum cutoff in angstrom
   * - cutoffmax
     - Maximum cutoff in angstrom
   * - symf_calc_or_NN
     - symf_calc if you want to calculate only ithe symmetry function. NN if you want to calculate the NN result
   * - NNrdf
     - True if you want to calculate the distance from surface
   * - fullNeigh
     - Number of neighbors within the cutoff. If molecule has fewer neighbors, it will be concidered as surface

   
The following is the parameter information about the symmetry function. This part functions as a block, when you edit the INPUT.dat file, make sure you have them in the right order. 
     
.. list-table::
   :widths: 25 70
   :header-rows: 1
  
   * - Name 
     - Description
   * - number_of_para
     - You need to have the number of a specific symmetry function you are using. 
   * - point
     - If you use point representation, set as TRUE. If you use point-vector representation, set as FALSE.
   * - G2G3
     - If you set TRUE in the "point" column, this should be the next line. Choose a type of symmetry function: G2 or G3.
   * - vec_id
     - If you choose FALSE in the "point" column, write which element id for the atom you use for the vector. These element IDs must be unique in the molecule.   
   * - Rs
     - If you choose G2 in "G2G3", write parameters for Rs here separated by space. 
   * - eta
     - If you choose G2 in "G2G3", write parameters for eta here separated by space. 
   * - kappa
     - If you choose G3 in "G2G3", write parameters for kappa here separated by space. 


Output files
""""""""""""

Depending on which mode you use, you will get the different output. 

.. list-table::
   :widths: 25 70
   :header-rows: 1
  
   * - file name
     - Description
   * - atominfo.txt
     - Print out molecule type, element type, and coordinate information. 
   * - symf_result.txt
     - If you slected symf_calc, this file provides the symmetry function output.
   * - new_traj.dump
     - Print out trajectory data with Neural Network output.
   * - run_ave_traj.dump
     - Print out running average of Neural Network output with trajectory data. 
   * - NNrdf_result.dump
     - In addition to run_ave_traj.dump, it contains the number of local neighbor and distance to the surface


ver2.0
^^^^^^

This is version 2.0 of the Symmetry function codes. In this version, we use a different atom style called full2. providing more flexibility with input molecules. You can still use the version 1.0 code here as well. Ensure you use the right file type. If you are using this version, make sure to also read `USER-ATOM-VEC-FULL2 <USER-ATOM-VEC-FULL2_>`_. 

When using this version, make sure you have "sym_type" in trajectory file as follows:

For version 1.0:

ITEM: ATOMS id mol element x y z vx vy vz

For version 2.0:

ITEM: ATOMS id mol element sym_type x y z vx vy vz

.. Note::

   If you want to create a dump file using LAMMPS, you can use the following LAMMPS extension to insert "sym_type": `USER-ATOM-VEC-FULL2 <USER-ATOM-VEC-FULL2_>`_.


You have two options for the new style. Keep in mind that you can still use the old version too.

.. code-block:: console

   ./sa_test trajectory_file file_type  name_of_input_data.data

List of available file types:

.. list-table:: 
   :widths: 25 70
   :header-rows: 1

   * - File Type
     - Description
   * - lammps_vec
     - Use this command if you have an orthogonal box
   * - lammps_vec_tri
     - Use this command If you have a non-orthogonal box
   * - lammps_sym_new
     - NEW: Use this command if you use new version with an orthogonal box
   * - lammps_sym_new_tri
     - NEW: Use this command if you use new version with a triclinic box




USER-ATOM-VEC-FULL2
^^^^^^^^^^^^^^^^^^^


Here we introduce a new atom style in LAMMPS called 'full2'. 'full2' is similar to the 'full' atom style, which stores molecular parameters (bonds, angles, dihedrals, and impropers) and charge. Additionally, 'full2' stores symmetry function IDs. If you wish to use this atom_style, you need to have the following data file in the specified order:

.. code-block:: console

     Atom_ID  Molecule_ID  Atom_type  Charge  Sym-type  pos image_data(optional)

Output
""""""

The 'dump_custom' command can output the dump file with symmetry function IDs. To enable this output, include the following line in the run file, where 'traj.dump' will be the name of the output file:

.. code-block:: console

   dump              d1 all custom 1000 traj.dump id mol element sym_type x y z 


Libraly code with LAMMPS
^^^^^^^^^^^^^^^^^^^^^^^^

We have created several files to calculate the symmetry functions. If your system is larger than the current setting, you should modify line 12 and 14 in dafed.h.

After modifying any of the files, compile the code by typing: 

.. code-block:: console

   make

To run the command, make sure to load the Intel and OpenMPI modules. If you encounter any error messages, there may be a problem in your code that needs to be fixed. 

If you connect the libraly code with the LAMMPS source code, you need to copy dafed.h to dafed_link.h and delete the first line of dafed_link.h.

After fixing all the bugs, proceed as follows:

.. code-block:: console

   make lib

to compile the new changes to the LAMMPS source code. If you've added a new file, make sure to add it to the Makefile. 

Lammps source code
------------------

To compile the newest version of the LAMMPS code, go to the src directory and type:


.. code-block:: console

   make mpi 

If you have a new file to add in the source directory, t is I highly recommended you to create a subdirectory within the source directory and copy the file there. 
If you encounter any error messages after typing ``make mpi``, carefully review and fix all the problems you encounter. 

Running LAMMPS code
^^^^^^^^^^^^^^^^^^^

To run the LAMMPS code, you need to have a run file and the input files required in the run file. For instance, you need to have weight.txt, an input file, urea.settings, and urea.init file for running example.

Use USER-G2G3
^^^^^^^^^^^^^

This code can calculate the symmetry function using LAMMPS. Currently, this code accepts up to 2 point vector representations at maximum. If you would like to use this feature, make sure to copy the .cpp file and .h file to the src directory in LAMMPS and compile it.
The following line is the syntax for the LAMMPS command: 

.. code-block:: console

   compute  ID  group-ID  G2G3gen/vecchunk center-atom-id  1st-vector-ID  2nd-vector-ID cutoff number_of_symmetry_function


Data Preparation
""""""""""""""""

In addition to the run file and input files, you also need INPUT.dat file. 

.. list-table::
   :widths: 25 70
   :header-rows: 1
  
   * - name 
     - Description
   * - G2 point / G3 point
     - Write number of G2 or G3 point symmetry function.  
   * - G2p Rs / G2p eta 
     - Write Rs or eta for the G2 point symmetry function.
   * - G3p kappa
     - Write kappa for the G3 point symmetry function.
   * - G2 point vec1 / vec2
     - Write number of G2 point vector symmetry function.
   * - G2v1 Rs / G2v2 Rs 
     - Write Rs for the G2 point vector symmetry function.
   * - G2v1 eta / G2v2 eta
     - Write eta for the G2 point vector symmetry function.   
   * - G3 point vec1 / vec2 
     - Write number of G2 point vector symmetry function.
   * - G3v1 kappa / G3v2 kappa
     - Write kappa for the G3 point vector symmetry function.


Use USER-NNout and USER-NNoutg2g3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This code can calculate the classification result using NN with symmetry function parameter. If you would like to use this, make sure you copy the .cpp file and .h file to src directory in LAMMPS and compile it.

Data Preparation
""""""""""""""""

This file requires an INPUT.dat file, which has a similar style to the Symmetry function codes, but with slight modifications. Additionally, you need to have weight.txt, which you can prepare in the `Neural Network Classification section <file:///home/drk354/Downloads/MolStrucClassifier-main/docs/build/html/neuralnet.html>`_. The "Noutg2g3" code can store misclassified molecules in addition to the symmetry function result. 

Basic parameters Required in INPUT.dat:

.. list-table::
   :widths: 25 70
   :header-rows: 1

   * - Name 
     - Description
   * - num_sf
     - Number of total symmetry functions
   * - num_pvsf_type
     - Number of point-vector symmetry functions
   * - COvectype
     - Element IDs for the first vector 
   * - NNvectype
     - Element IDs for the second vector 
   * - triclinic
     - True if this is in a non-orthogonal box
   * - center atom
     - Element ID of the center atom. Must be unique in the molecule
   * - atm in molec
     - Number of atoms in the molecule. If you have ia defect, choose the one with more atoms 
   * - cutoffmin
     - Minimum cutoff in angstrom
   * - cutoffmax
     - Maximum cutoff in angstrom

Parameter information about the symmetry function is the same as the `Parameters section of ver1.0 <file:///home/drk354/Downloads/MolStrucClassifier-main/docs/build/html/How_to_use.html#parameters>`_
In the run file, you need to specify the element ID of the center atom and the number of  symmetry function. "NNoutg2g3" requires an additional argument called "prob_error". If "prob_error" is 0.4, it means that if less than 40% of the local environment is classified the same as the targeted molecule. 

Syntax for NNout is:

.. code-block:: console

   compute CNN all nnoutg2g3/atom center_id #of_symmetry_function

Syntax for NNoutg2g3 is:

.. code-block:: console

   compute CNN  all nnoutg2g3gen/atom center_id #of_symmetry_function prob_error

Output files
""""""""""""

The output of "NNout" is stored in "g2g3output.dat", which contains both the symmetry function and classification results. 
If you use "NNoutg2g3" instead, additional files will be generated. "g2g3good.dat" stores the symmetry function and classification results for correctly classified molecules, while "g2g3error.dat" contains similar information for molecules that are misclassified.





How to use Ovito
----------------

I used Ovito for visualizing trajectory files and add some necessary information such as mass of atoms etc.

Step 1
^^^^^^

Assume you were able to download Ovito. If not, go to `Ovito webpage <https://www.ovito.org/windows-downloads/>`_ and download it. You need the free version for this project. 

Step 2
^^^^^^

In Ovito, you have two ways to read files. You can read files from your computer or read files remotely. When you want to use SSH to connect, make sure you add ``sftp://`` at the beginning. After you read the file, you can visualize the file and work on it. If you see error messages, you may have wrong data type/file or lost connection to the remote directory. In the latter case, close the window and open Ovito again. 

Use Ovito with VASP file
------------------------

We will focusing on VASP file in this page, but you can use other data file as long as the final product has the same atom style as LAMMPS' full style. After you get the VASP file from Ogre, you can use Ovito to make sure the data did not corrupt. If you want to make an interface, you can check 
:ref:`How to use Ogre <How_to_use_Ogre>`

In this page, I will focus on what you should do in general to make a LAMMPS data file. 

If you are working on a urea molecule, you need to add the following information.

Add Masses
^^^^^^^^^^

First, you need to add particle types from Data Source on the right part. You can choose whatever name you want for the name, but I usually use hn, o, c, n1, and n2. The masses of each atom are c = 12.01, hn = 1.008, n = 14.01, o = 16.00. 

Add Charge
^^^^^^^^^^^

Select Type in Add Modification.

choose the atom type.

Compute property in Add Modification.

Change output property to Charge.

Click `Compute only for selected elements`.

Choose values for Expression (c = 0.142, hn = 0.333, n = -0.542, o = -0.39).


Modify Positions
^^^^^^^^^^^^^^^^

1. Click Manual Selection From Add modification

2. Pick atoms that you want to modify, if you select something before, don't forget to clear the selection

3. Click Compute Property from Add modification

4. Change output property to Position

5. Click Compute only for selected elements

6. Initially, choose Position.X,Y,Z

7. Then add or subtract CellSize.X,Y,Z

8. Repeat this until all the position get fixed

9. If you are working on interface, basic move is same the as above, BUT when you modify atoms from the different cell size, you have to use original Cell Size to shift the molecule.


Add Bonds, cluster, and molecule identifier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Choose Create bonds for all types using a pairwise cutoff (typically 1.4).

2. Select Cluster analysis from Add modification.

3. Choose Bonds.

4. Click Compute property in Add modification.

5. Type Molecule identifier for output property and Cluster for Expression.

6. Go to Compute Property again and type Periodic image. Choose 0, 0, 0 for xyz.

7. Export as data file or VASP file.


How to Use moltemplate
----------------------

Below are the basic steps for working with Moltemplate. I'll explain one of the examples below:

If you haven't worked with Moltemplate recently, get system.lt and urea.lt to the test directory.

Run the following command for lammps run:

.. code-block:: console
   
   moltemplate.sh system.lt

You can get rid of unnecessary lines in the file by typing the following command:

.. code-block:: console
   
   cleanup_moltemplate.sh

You will get these files after running the previous commands:



.. table:: 
   :align: left

   +----------------+----------------------+
   | system.in      |lammps input file     |
   +----------------+----------------------+
   | system.data    |lammps data file      |
   +----------------+----------------------+
   | system.in.init |force field functions |
   +----------------+----------------------+

If you finish the urea sample, you can continue with the actual data. 

Step 1
^^^^^^

Get the file from ovito with mass, charges, and bond information (in the format of a .data file). If you have other information, you can delete it now. Then type:

.. code-block:: console
   
   ltemplyfy.sh ovitofilename > outputfilename     (ltemplyfy.sh urea.data > urea.lt)

For example:

.. code-block:: console

   ltemplyfy.sh urea.data > urea.lt
   


Step 2
^^^^^^

Check the urea.lt file or the output file you made. You may have to modify several things. If you don't see a part that I will mention, skip that part.

1. Add import "gaff.lt" before write_once("In Init") {.

2. Add Urea inherits GAFF { on the next line after the import "gaff.lt".

3. Add another closed curly bracket } at the end of the file.

4. If you have a unique name for the atom type that is not used in GAFF, make sure you change the name in the Data Atoms section. For example, if you have n1 and n2, you can keep the name in Data Masses, but you have to have the right name such as n in Data atom section. Check useful commands shown below.

5. Remove bond type (such as @bond:type1).

6. Change "Data Bonds" to "Data Bond List".  

7. Check the atom information. If it says type##, you have to change the name of the type to the corresponding one. If you change this, you also have to change @atom: to the new name. Also you need to change atom in Data Bond List.

8. In Data Atoms section, if you don't have $mol section or the value for this is the same for the whole system, you have to check the input data. Most likely you need to go to Ovito and assign molecule data. Basically, you need to do clustering and go to compute property > molecule identifier > cluster. More details are in How to use Ovito.

9. If you see angle data or some others after Bond information, remove that information now.

10. module load anaconda3/2020.07

11. Change system.lt's first line to whatever the input file, change the cell box size.

12. Run:

.. code-block:: 

   moltemplate.sh system.lt

13. Check system.data to make sure it ran.

14. Clean up the code by typing:

.. code-block:: 

   cleanup_moltemplate.sh


15. If you want to change the atom type, do it now. Make sure you change it to the right name.

16. Change the n2 in the Atom List to the corresponding name.

17. Use this file as a standard.

18. When you use it in LAMMPS, you need to have system.data, system.in.init, and system.in.settings in the running directory.

19. You need to modify urea.in.settings if you change the atom information. For instance, if you added n2, you should have a line for that too.

20. Also, change the atom type from 4 to 5 for urea in system.data.


Useful Commands When Using Moltemplate in Vim:


.. code-block:: 

   :/@bond:type3           <-- highlight the words type3.

.. code-block:: 

   :%s///g                 <-- delete the highlighted words.

.. code-block::

   %s/@atom:n1 /@atom:n /g	<-- replace all @atom:n1 in the file to @atom:n
   
   
.. _How_to_use_Ogre:


How to use Ogre
-----------------------

If you have any questions about how to install or fix the initial error, check the install ``Ogre`` page for the detailed instructions. 

Step 1
^^^^^^

For all the parameters, you can check the ogre.config file. This is the only file you need to change if the program is running. The following is the sample config setting for my research: 

.. code-block:: console

   [io] 
   structure_path = urea/urea_poly4.vasp
   ; THis is where the original file is stored at. I used a vasp file which was generated from ogre. 
   structure_name = UREA_1_4_for_interface
   ;The new files will be stored in this directory
   format = VASP
   ; Format can be FHI, VASP or CIF. I used VASP file for the initial file. 
   [methods]
   cleave_option = 0
   ;0: cleave a single surface, 1: cleave surfaces for surface energy calculations
   [parameters]
   layers = 1-5
   ; This tells the number of the layers you would like to have as an output. In this case, you will get 1 layer to 5 layers in the separate file. 
   vacuum_size = 40
   ;This tells the size of the vacuum region in z access. The box size in z direction will be 40+40 + actual box size
   highest_index = 3
   ; Only needed in cleave_option = 1， ignored in cleave_option = 0
   supercell_size = None
   miller_index = 1 0 0
   ; Only needed in cleave_option = 0, ignored in cleave_option = 1. This changes the direction of where you get the vacuum layer
   desired_num_of_molecules_oneLayer = 0
   ; Set to 0 if you don't want to change one layer structure. Default is 0 (highly recommended)
   
   
After you finish the parameter set up, you can run the program by typing

.. code-block:: console

   python runOgre.py --filename ogre.config

.. warning::

   When you use the vasp file from ogre and modify the data in ovito, make sure you delete the first line of the vasp file. This data is redundant and Ovito will crash if it is there. 
   You can use the POSCAR file as an output file. If you want to make an interface of two different poscar files, you need to do the following. 

Step 1
^^^^^^^

You need to use atoms_shifter.py to shift a file. Change your input file to the appropriate name and change the amount you need to shift. 

Step 2
^^^^^^

Open up the shifted data and another poscar file you wish to combine. Copy line 6 and 7 and add them right next to the other file's line 6 and 7. DO NOT add numbers on line 7 this will mess all the data up. 
THis is one example of an interface I made. First half was from Urea Poly I and the second half came from Urea Poly IV. 
H C S N O H C N S O
1920 480 480 480 480 1440 360 360 360 360

After adding line 6 and 7, add the second interface atoms coordinate data after the first one. 
















