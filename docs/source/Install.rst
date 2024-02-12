.. include:: .colors.rst

Install
=======

This page will explain how to install all the packages you need for this project step by step. 

We use NYU high performance computer GREENE as an example. If you do not have access to the account, go to `NYU HPC account page <https://sites.google.com/nyu.edu/nyu-hpc/accessing-hpc/getting-and-renewing-an-account?authuser=0>`_ to get an access. 

.. _installation:




Install LAMMPS libraly code
---------------------------

Step 1
^^^^^^

First of all, you need a github data from Jutta Rogal's github repository. If you do not have an access to it, contact ''jr4855@nyu.edu'' to get the access.

After you get an access, type the following command. ``This step is only for the first time user``

.. code:: console

   git clone https://github.com/rogalj/enhanced_sampling 

If you already have an access and installed the data, you can update the information by 

.. code:: console

   git pull

If you want to check if the file is up to date or not, you can check by typing 

.. code:: console

   git update

If you want to add some directories/files to the gthub repository, you can upload them by typing

.. code:: console

   git add name_of_folder/file

This process is not in repository yet.  ``git status`` is one of the useful commands that can check the status of the local repository. This command will show all the files that you added, deleted or modified which are not yet reflected to the local repository. 

When you upload the directory/files, you can also make an comment. You can add them like the following

.. code:: console

   git commit -m "WHATEVER THE COMMENT YOU WANT"

This process will add to your local repository. If you want to add them to the global repository, type 

.. code:: console

   git push

At this point, your work will be on github repository.

IF you made a mistake and added unnecessary file,

.. code:: console

   git reset

to cancel the command you were working on.




Install Lammps code
-------------------

In order to run the libraly programs, you need to have LAMMPS code. You can search how to get from internet or follow this.

Step 1
^^^^^^

Go to the directory where you want to store the lammps code. I use mylammps as an example. Run the following command at mylammps directory to install Lammps from github.

Run the following command

.. code:: console

    git clone -b release https://github.com/lammps/lammps.git mylammps

If you want to add files (such as fix, commpute files), add them in the src directory. If you want to add directory, you can add them at source directory as well. 


Step 2
^^^^^^

Next step is getting an appropriate module into your system. If you are not working on NYU hpc, search a newest version of intel and openmpi. 

.. code:: console

   module load intel/19.1.2
   module load openmpi/intel/4.1.1

For intel, intel/17.0.0.098 or later version should work. 

For openmpi, openmpi/2.1.6 or later version should work.  

After you load, update any information in src directory, make sure you update them by typing 

.. code:: console

   make mpi

Step 3
^^^^^^

The next step is run the lammps code with using Symbolic link and soft link. (wheich enable to run the program from any of your directory) 


On GREENE, people usually make a link in a bin directory. Let's make a bin directory if you do not have one. (/home/`username`/bin). 

.. code:: console

   ln -s /home/:green:`username`/mylammps/src/lmp_mpi lammps

Make another softlink that connect the lammps code and the libraly files. Make a file in libraly directory like Dafed_DK.inc In this file you should have two things to connect the file. 
When you make this file, make sure you cange the path to the appropriate name for your case. 
   

.. code:: console

   DAFED_LOAD= /home/:green:`username`/ES_result/molvec/libcv_nn.so -ldl
   DAFED_DEPENDENCIES = /home/:green:'username'/ES_result/molvec/libcv_nn.so
   
After this, go to mylammps directory and type the following commands 

.. code:: console

   ln -s /home/:green:'username'/ES_result/molvec/Dafed_DK.inc Dafed.inc at the parent directory of lammps src
   ln -s /home/:green:'username'/ES_result/molvec/dafed_link.h Dafed.h at the parent directory of lammps src


Step 4
^^^^^^

You are almost there! The next step is install some additional lammps package at the source directory. You will need the following packages. 

.. table:: 
   :align: left

   +---------------+
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

The command you need to type to get these additional packages are the folowing

.. code:: console

   make yes-:green:'name_of_package'

You can check the status from typing 

.. code:: console

   make pi

Open ``Makefile.package`` and add the following to connect the source directory to the libraly directories. 

.. code:: console

   PKG_LIB =    $(DAFED_LOAD)



Install moltemplate
-------------------

In our programs, you need to have specific data file called .data This is a common data type for lammps programs, however, you may not be able to find this data set online. Other common dataset are .vasp fole, .cif, .xyz file etc. If you have H.data type with missing mose critical information (such as bonds, angles, dihedrals, improper angles infomation), you can use moltemplate to fill these. If you do not have .data file, check :doc:`How_to_use` for more information. 

Step 1
^^^^^^

run the following command to install moltemplate

.. code:: console

   git clone https://github.com/jewettaij/moltemplate moltemplate


Step 2
^^^^^^

To create a symbolic link, you need to go to bin file (/home/`username`/bin) as shown in Step3 of Installing Lammps. 

.. code:: console

   git clone https://github.com/jewettaij/moltemplate moltemplate

Step 3
^^^^^^

Python3.7.2 or above is required to use for this program. In greene, we have python 3.8.6 and you have to load this module like the following.

.. code:: console

   module load python/intel/3.8.6

Step 4
^^^^^^

go to moltemplate directory and finish installation

.. code:: console
   
   pip3 install .



Install Ogre
------------

Ogre is a computational tool to create a surface with vacuum region. In project, I used this code to create surfaces. 


Step 1
^^^^^^

First of all, if you want to download the code you need to go to `Ogre webpage <https://www.noamarom.com/software/download/>`_ and choose Ogre2.0.


Step 2
^^^^^^

After you download the program, there are some bugs in the code that you have to fix. Go to ``/Ogre-master/ogre/utils`` and remove the following two lines. 

.. code:: console
   
   from ogre.utils.surface_energy import convergence_plots
   from ogre.utils.unique_planes import UniquePlanes

Go to ``/Ogre-master/ogre/utils`` and remove the following two lines. 

.. code:: console
   
   from ogre.utils.surface_energy import convergence_plots
   from ogre.utils.unique_planes import UniquePlanes

Go to ``/Ogre-master/ibslib/io`` and open ``check.py``. 

change the line 22 from 


.. code:: console
   
   extension2format = ase_extension2format`
  
to 

.. code:: console

   extension2format = deepcopy(ase_extension2format)









