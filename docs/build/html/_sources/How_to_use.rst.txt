How to link with LAMMPS and run
===============================


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
^^^^^^^^^^^^^^^^^^^^

To run the LAMMPS code, you need to have a run file and the input files required in the run file. For instance, you need to have weight.txt, an input file, urea.settings, and urea.init file for running example.

Symmetry function calculation with LAMMPS (USER-G2G3)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


Neural Network calculation with LAMMPS (Use USER-NNout and USER-NNoutg2g3)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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












