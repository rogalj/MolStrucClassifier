.. include:: .colors.rst

How To Install
===============

This page will explain how to install all the packages you need for this project step by step. 



.. _installation:

Requirement
-----------

This package requires C++ compilar. Environments under intel 19.1.2 was tested


Install libraly code (MolStrucClassifier)
--------------------------------------------

Copy the repository
^^^^^^^^^^^^^^^^^^^

You can copy the repository by typing the following command:

.. code:: console

   git clone https://github.com/rogalj/MolStrucClassifier.git

If you already installed the repository, you can update it by typing:

.. code:: console

   git pull

If you are an editter, check the following command as well. 

If you make any changes and want to add directories/files to the GitHub repository, you can upload them by typing:

.. code:: console

   git add name_of_folder/file

The following command can check the status of the local repository. This command will show all the files that you added, deleted or modified which are not yet reflected to the local repository. 

.. code:: console

   git status
   

When you upload the directories/files, you can also add a comment. You can add them like the following:

.. code:: console

   git commit -m "WHATEVER THE COMMENT YOU WANT"

This process will add to your local repository. If you want to add them to the global repository, type:

.. code:: console

   git push

At this point, your work will be on the GitHub repository.


If you made a mistake and added an unnecessary file, type:


.. code:: console

   git reset

to cancel the command you were working on.




Install LAMMPS code
-------------------

In order to run the libraly programs in LAMMPS, you need to have the LAMMPS code. If you don't need to connect with LAMMPS, you can skip this section. If you already have LAMMPS, make sure you connect the libraly and have the appropriate packages shown in Step 4. 

Step 1
^^^^^^

Go to the directory where you want to store the lammps code. I use mylammps as an example. Run the following command at mylammps directory to install LAMMPS from github.

Run the following command

.. code:: console

    git clone -b release https://github.com/lammps/lammps.git mylammps

If you want to add files (such as fix, compute files), add them in the src directory. You can also add directories in the src directory.




Step 2
^^^^^^

The Next step is to get an appropriate module into your system: 

.. code:: console

   module load intel/19.1.2
   module load openmpi/intel/4.1.1
   module load gcc/10.2.0
   
We also tested with Intel 17.0.0.098 and OpenMPI 2.1.6.  

After you load modules, update files in the src directory by typing:


.. code:: console

   make mpi

Step 3
^^^^^^

The next step is to run the LAMMPS code using symbolic links and soft links. (which enable you to run the program from different directories) 


On GREENE, people usually make a link in a bin directory. Let's make a bin directory if you don't have one. (first go to /home/username/bin). 

.. code:: console

    ln -s /home/username/mylammps/src/lmp_mpi lammps

Make another softlink that connects the LAMMPS code and the libraly files. Make a file in the libraly directory like `Dafed_DK.inc`. In this file you should have two things to connect the file. 
When you make this file, make sure to change the path to the appropriate name for your case. 
   

.. code:: console

   DAFED_LOAD= /home/drk354/MolStrucClassifier/MolStrucClassifier/molvec_lib/ver2.0/libcv_nn.so -ldl
   DAFED_DEPENDENCIES = /home/drk354/MolStrucClassifier/MolStrucClassifier/molvec_lib/ver2.0/libcv_nn.so
   
After this, go to the parent directory and type the following commands:

.. code:: console

   ln -s /home/drk354/MolStrucClassifier/MolStrucClassifier/molvec_lib/ver2.0/Dafed_DK.inc Dafed.inc 
   ln -s /home/drk354/MolStrucClassifier/MolStrucClassifier/molvec_lib/ver2.0/dafed_link.h Dafed.h 


Step 4
^^^^^^

You are almost there! The next step is to install some additional LAMMPS packages at the source directory. You will need the following packages:

.. table:: 
   :align: left

   +---------------+
   |name_of_package|
   +===============+
   | KSPACE        |
   +---------------+
   | MANYBODY      |
   +---------------+
   | MOLECULE      |
   +---------------+
   | OPT           |
   +---------------+
   |EXTRA-MOLECULE |
   +---------------+

The command you need to type to get these additional packages is as folows:

.. code:: console

   make yes-name_of_package

You can check the status from typing:

.. code:: console

   make pi

Open ``Makefile.package`` and add the following to connect the source directory to the libraly directories:

.. code:: console

   PKG_LIB =    $(DAFED_LOAD)

Open ``Makefile.package.setting`` and add the following to connect at the end of the file:

.. code:: console

   include ../../Dafed.inc



Install moltemplate
-------------------

In order to run the libraly programs in LAMMPS, you need to have specific data file called .data. If you cannot find the dataset online, check the following section to create one. If you have a .data file with missing critical information (such as bonds, angles, dihedrals, improper angles information), you can use moltemplate to fill these. If you do not have .data file, refer to :doc:`How_to_use` for more information. 

Step 1
^^^^^^

Run the following command to install Moltemplate:

.. code:: console

   git clone https://github.com/jewettaij/moltemplate moltemplate


Step 2
^^^^^^

To create a symbolic link, navigate to the bin directory (/home/username/bin), as shown in Step 3 of Installing LAMMPS. 

.. code:: console

   git clone https://github.com/jewettaij/moltemplate moltemplate

Step 3
^^^^^^

Python3.7.2 or above is required to use for this program. On Greene, we have python 3.8.6 and you have to load this module as follows:

.. code:: console

   module load python/intel/3.8.6

Step 4
^^^^^^

Navigate to the Moltemplate directory and finish the installation:

.. code:: console
   
   pip3 install .



Install Ogre
------------

Ogre is a computational tool to create a surface with a vacuum region. In this project, I used this code to create surfaces. 

Step 1
^^^^^^

If you want to download the code, go to `Ogre webpage <https://www.noamarom.com/software/download/>`_ and choose Ogre2.0.

Step 2
^^^^^^

After downloading the program, make some minor changes if the code does not work. Go to ``/Ogre-master/ogre/utils`` and remove the following two lines:


.. code:: console
   
   from ogre.utils.surface_energy import convergence_plots
   from ogre.utils.unique_planes import UniquePlanes

Next, go to ``/Ogre-master/ogre/utils`` and remove the following two lines:

.. code:: console
   
   from ogre.utils.surface_energy import convergence_plots
   from ogre.utils.unique_planes import UniquePlanes

Go to ``/Ogre-master/ibslib/io`` and open ``check.py`` and change line 22 from:

.. code:: console
   
   extension2format = ase_extension2format`
  
to:

.. code:: console

   extension2format = deepcopy(ase_extension2format)


Using NYU GREENE
-----------------

We will use NYU's high-performance computer, GREENE, as an example. If you do not have access to the account, visit the `NYU HPC account page <https://sites.google.com/nyu.edu/nyu-hpc/accessing-hpc/getting-and-renewing-an-account?authuser=0>`_  to request access or download the appropriate softwares on your computer.






