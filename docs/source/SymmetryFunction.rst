.. _SymmetryFunction:

Symmetry Function codes (molvec_lib)
====================================

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

This is version 2.0 of the Symmetry function codes. In this version, we use a different atom style called full2. providing more flexibility with input molecules. You can still use the version 1.0 code here as well. Ensure you use the right file type. If you are using this version, make sure to also read :ref:`FULL2`.

Data Preparation
""""""""""""""""

When using this version, make sure you have "sym_type" in trajectory file as follows:

For version 1.0:

ITEM: ATOMS id mol element x y z vx vy vz

For version 2.0:

ITEM: ATOMS id mol element sym_type x y z vx vy vz

.. Note::

   If you want to create a dump file using LAMMPS, you can use the following LAMMPS extension to insert "sym_type"\: :ref:`FULL2`.


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

Output files
""""""""""""""""

You will get the same output as shown in ver1.0. 
