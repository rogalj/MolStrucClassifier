.. _NYU:

For NYU HPC users
==================






Symmetry Function and LAMMPS code
""""""""""""""""""""""""""""""""""

We will use NYU's high-performance computer, GREENE, as an example. If you do not have access to the account, visit the `NYU HPC account page <https://sites.google.com/nyu.edu/nyu-hpc/accessing-hpc/getting-and-renewing-an-account?authuser=0>`_  to request access or download the appropriate softwares on your computer.


Greene computing cluster
++++++++++++++++++++++++++++

On Greene, we tested with the following modules. 

.. code:: console

	module load intel/19.1.2
	module load openmpi/intel/4.1.1

After 2024 March HPC update, we also need to use GCC module as well. If you cannot compile symmery function files, make sure you load the following module as well. 

.. code:: console

	module load gcc/10.2.0


Rusalka computing cluster
+++++++++++++++++++++++++++

On Rusalka, we tested on older version of intel and openmpi. If you load the following modules, it should work. 

.. code:: console

    module load intel/17.0.0.098
    module load openmpi/3.1.3
    
    
    
Neural Network Classification Code
""""""""""""""""""""""""""""""""""

Greene computing cluster
+++++++++++++++++++++++++++

On Greene, all the necessary Intel MKL libraries and compilers should be installed.  They are accessible by loading the corresponding modules, e.g.

.. code:: console

    module load openmpi/intel/4.1.1

Then the provided ``Makefile`` can be use and the code is compiled by

.. code:: console

   make

which will produce an executable ``ml-classification.exe``.  This executable will run in parallel with the number of threads given in the input file.


Rusalka computing cluster
+++++++++++++++++++++++++++

On Rusalka, all the necessary Intel MKL libraries and compilers should be installed.  They are accessible by loading the corresponding modules, e.g.

.. code:: console

    module load intel/17.0.0.098
    module load openmpi/3.1.3

Then the provided ``Makefile`` can be use and the code is compiled by

.. code:: console

   make

which will produce an executable ``ml-classification.exe``.  This executable will run in parallel with the number of threads given in the input file.

