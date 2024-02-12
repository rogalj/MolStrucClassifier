How to use
==========



Lammps libraly code
-------------------

I created main.cpp, molvec.cpp, and dafed.h. If you have a system that is larger than current setting, you should modify dafed.h line 12 and 14.

After you modify any of the file, you need to update the information to share via typing 

.. code-block:: console

   make

In order to run the command, you need to load intel and openmpi modules. If you see any error message, you have a problem in your code and it won't work unless you fix it. 

When you connect with lammps source code, you need to copy dafed.h to dafed_link.h and delete the first line of dafed_link.h.

After you fix all the bugs you need to type 

.. code-block:: console

   make lib

to upload the new change to lammps source code. If you made a new file, make sure you add them to Makefile. In this file, you have to addt eh new file to all section and lib section. 

Lammps source code
------------------

When you want t oupload the newest version of lammps code, you need to type 

.. code-block:: console

   make mpi 

at the src directory. If have a new file to add in the source directory, I highly recommend you to make a subdirectory in the source file and copy that file to source directory. 
If you see any error message after typing ``make mpi``, take a closer check and fix all the problems you have. 

If you want to modify the code/parameters, change fix_dafed.cpp line 504-517, fix_dafed.h line 135, 136. If you are working with compute_NNout, change line 149-171 in cpp file, and change line 77,78 in h file. 



Running LAMMPS code
-------------------

In order to run the lammps code, you need to have run file, slurm file, any other input file required in the run file. For instance, you need to have weight.txt, input file, urea.settings, urea.init file for running my code. 


How to use Ovito
----------------

I used Ovito for visualizing trajectory files and add some necessary information such as mass of atoms etc.

Step 1
^^^^^^

Assume you were able to download Ovito. If not, go to `Ovito webpage <https://www.ovito.org/windows-downloads/>`_ and download it. You need a free version for this project. 

Step 2
^^^^^^

In Ovito, you have two ways to read files. You can read files from your computer or read file remotely. When you want to use ssh to connect, make sure you add ``sftp://`` at the beginning. After you read the file, you can visualize the file and able to work on it. If you see the error messages, you may have wrong data type/file or lost connection to the remote directory. In latter case, close the window and open Ovito again. 

Use Ovito with VASP file
------------------------

After you get the vasp file from Ogre, you can use Ovito to make sure the data did not corrupt. If you want to make an interface, you can check 
:ref:`How to use Ogre <How_to_use_Ogre>`

In this page, I will focus on what you should do in general to make a lammps data file. 

If you are working on urea molecule you need to add the following information.

Add masses
^^^^^^^^^^

First, you need to add particle type from Data Source on a right part. You can choose whatever the name you want for the name, but I usually use hn, o, c, n1, and n2. The masses of each atom are c = 12.01, hn = 1.008, n = 14.01, o = 16.00. 

Add Charge
^^^^^^^^^^^

Select type in Add modification

choose the atom type

Compute property in Add modification

change output property to Charge

click Compute only for selected elements

Choose values for Expression (c = 0.142, hn = 0.333, n = -0.542, o = -0.39)


Modify Positions
^^^^^^^^^^^^^^^^

Click Manual Selection From Add modification

Pick atoms that you want to modify, if you select something before, don't forget to clear the selection

Click Compute Property from Add modification

Change output property to Position

Click Compute only for selected elements

Initially, choose Position.X,Y,Z

Then add or subtract CellSize.X,Y,Z

Repeat this until all the position get fixed

If you are working on interface, basic move is same the as above, BUT when you modify atoms from the different cell size, you have to use original Cell Size to shift the molecule.


Add Bonds, cluster, and molecule identifier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Choose create bonds for all types using pair wise cutoff (typically 1.4)

Select Cluster analysis from Add modification

Choose Bonds

Click Compute property in Add modification

Type Molecule identifier to output property and type Cluster for Expression

Go to compute property again and type periodic image. choose 0, 0, 0 for xyz

Export as data file or vasp file.






How to Use moltemplate
----------------------

The following are the basic move on moltemplate. I will explain one of the example below, 

If you have not worked with moltemplate recently, get system.lt and urea.lt to the test directory.

Run the following command for lammps run.

.. code-block:: console
   
   moltemplate.sh system.lt

You can get rid of unnecessary lines in the file by typing the following command.
.. code-block:: console
   
   cleanup_moltemplate.sh

You will get these file after running the previous commands.

.. table:: 
   :align: left

   +---------------------------------------+
   | system.in (lammps input file)         |
   +---------------------------------------+
   | system.data (lammps data file)        |
   +---------------------------------------+
   | system.in.init (force field functions)|
   +---------------------------------------+

If you finish the urea sample, you can continue with the actual data. 

Step 1
^^^^^^

Get file from ovito with mass, charges, and bond information. (format of .data file) If you have other information, you can delete them now. Then type 

.. code-block:: console
   
   ltemplyfy.sh ovitofilename > outputfilename     (ltemplyfy.sh urea.data > urea.lt)

to get the ltemplyfy file. 


Step 2
^^^^^^

Check urea.lt file or the output file you made. You may have to modify several things. If you don't see a part that I will mention, skip that part.

1. Add   import "gaff.lt" before write_once("In Init") {

2. Add Urea inherits GAFF {    next line to the add import "gaff.lt"

3. Add another closed the curly bracket } at the end of the file

4. If you have a unique name for the atom type that is not used in GAFF, make sure you change the name in Data Atoms section. For example, if you have n1 and n2, you can keep the name in Data masses, but you have to have a right name such as n in Data atom section.Check useful commands shon below

5. remove bond type (such as @bond:type1)

6. Change "Data Bonds" to "Data Bond List"
  
7. Check the atom information. If it says type##, you have to change the name of the type to coresponding one. If you change this, you also have to change @atom: to the new name. Also you need to change atom in Data Bond List.
   
8. In Data Atoms section, if you don't have $mol section or the value for this is same for the whole system, you have to check the input data. Most likely you need to go to Ovito and assign molecule data. Basially you need to do clustering and go to compute property > molecule identifier > cluster. More details are in How to use Ovito

9.  If you see angle data or some others after Bond information, remove that information now.
   

10. module load anaconda/2019........

11. change system.lt's first line to whatever the input file, change the cell box size.

12. 

.. code-block:: 

   moltemplate.sh system.lt

13. check system.data to make sure it ran

14. Clean up the code by typing 

.. code-block:: 

   cleanup_moltemplate.sh


15. If you want to change the atom type, do it now. Make sure you change to the right name.

16. Change the n2 in the Atom List to corresponding name

17. Use this file as a standard.

18. When you use in lammps, you need to have system.data, system.in.init, and system.in.settings. to the running directory.

19. You need to modify urea.in.settings if you change the atom information. For instance, if you added n2, you should have a line for that too.

20. Also change atom type from 4 to 5 for urea in system.data.



Useful commands when you use moltemplate in vim

.. code-block:: 

   :/@bond:type3           <-- highlight the words type3.

.. code-block:: 

   :%s///g                 <-- delete the highlighted words.

.. code-block::

   %s/@atom:n1 /@atom:n /g	<-- replace all @atom:n1 in the file to @atom:n





.. _How_to_use_Ogre:


How to use Ogre
-----------------------


If you have any questions about how to install or how to fix the initial error, check the install ``Ogre`` page for the detailed instructions. 

Step 1
^^^^^^

For all the parameters, you can check ogre.config file. This is the only file you need to change if the program is running. The following is the sample config setting for my research. 

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
















