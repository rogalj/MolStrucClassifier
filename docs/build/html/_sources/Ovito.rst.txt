.. _Ovito:

Get data file using ovito 
==========================


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

