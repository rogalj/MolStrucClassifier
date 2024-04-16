.. _neuralnet:

Neural Network Classification Code
==================================

The neural network classification code uses a feed-forward neural network to assign an output vector with a given value for each class it has been trained for.  The code can be used to do both training and testing.

The code can be found in the folder

.. code:: console

    MolStrucClassifier/nn-fit



The code
^^^^^^^^^^^^^^^^^^^^

The source code is ``nn-class.cpp`` together with the header file ``nn-class.h`` and can be compiled using the ``Makefile``, for details see `Compilation`_.  The code runs in parallel using mpi and depending on the platform and compiler the corresponding libraries need to be available accordingly.  The number of mpi threads is set in the input file, see description of the `input file <Main input file_>`_.

Within the folder, there is also an example for the input files ``input.dat`` and a training data set ``train.dat``, as well as an example file for a submission script to the `slurm queueing system <https://slurm.schedmd.com/overview.html>`_ ``job_NN.slurm``.


Hard-coded variables
"""""""""""""""""""""

The hard-coded variables can be found in the header file ``nn-class.h``.  The hard-coded variables include:

:MAX_NODES:  largest value of a hidden layer dimension + 1 (default is 51)

Furthermore, the following variable should be adjusted as needed:

:goodrange[2]{1,1}: scaling range for input and output values; this is an absolute value; the absolute value of the range has to be the same for input and output 

There is currently one activation functions implemented, 1/(1+exp(x)); again this is hard-coded and other activation functions can be added as needed. 



Input and output files 
^^^^^^^^^^^^^^^^^^^^^^^^

Main input file
"""""""""""""""""

The main input file contains the basic set-up; an example file ``input.dat`` is provided together with the code and shown in the following:

.. code:: console 

        ------------------------------------------------------------------------------------
         Input File for neural network fit
        ------------------------------------------------------------------------------------
        # initial seed for random number generator (>0)
        # if seed <= 0 a random seed between 1 and 1000 is generated
        ran_seed          123	

        # network architecture
        hidden layers     2 		! number of hidden layers
        input             14 		! number of input features
        hidden nodes      25 25  	! number of nodes in each hidden layer (one number for each layer)
        output            5 		! number of outputs

        # parallel threads
        nthreads          8

        # simulation setup
        run option        0        	! 1=train with new weights, -1=train with old weight; 0=test job
        iterations        250         	! maximum number of training iterations
        optimizer max     200    	! maximum number of optimizer steps
        save frequency    1  		! how often the results are saved

        train input       train.dat     ! name of file containing training data
        test input        test_A15.dat      ! name of file containing test data
        test output       output.dat    ! name of file for test output
        test error        error.dat     ! name of file for test error result (misclassified)


        weights input     weights_old.dat   ! name of file containing previous weights
        weights output    weights.dat   ! name of file to print new fitted weights
        cost output       cost.dat      ! name of the file to save values of cost function file fitting

        !-----------------------------------------------------------------------------------
        ! FORMAT OF THE INPUT FILE:
        ! the first 3 lines are simply read in as string
        ! the rest is parsed and identified over a string
        ! !IMPORTANT!  The first 16 (!) characters are reserved for the key word!!
        ! after this different parameter in one line must be separated by at least on space
        ! lines starting with '!' or '#' are comments




The first 3 lines are simply read in as string and ignored.
The first 16(!) characters of each line are reserved for the keyword and should not be changed; the values given for each keyword are free format; if several values are required, they need to be separated by spaces.  The order of the input lines can be changed.


:cost output:
   name of the output file containing the values of the cost function during the fit; (default: cost.dat)

:hidden layers:
   Number of hidden layers in the NN; this needs to be specified before the other parameters of the network architecture; (default: -1)

:hidden nodes:
   number of nodes in each hidden layer; provide one number of each layer separated by spaces;

:input:
   number of input features; (default: -1)

:iterations:
   number of optimisation steps; each step goes over all training points; (default: 1)

:nthreads:
   number of OMP parallel threads; (default: 1)

:optimizer max:
   maximum number of optimizer/gradient steps in each iteration; (default: 1)

:output:
   number of output nodes; (default: -1)

:ran_seed:
   Seed to initialize the random number generator; random numbers are used to initialize the weights of the NN in fitting; if a value smaller than 0 is chosen, a random seed between 1 and 1000 is set by the programme; (default: 42)

:run option:
   determines simulation task: 0 - read in weights file and evaluate test data; 1 - train NN with randomly initialized weights; -1 - train NN with weights read in from weights file; (default: 0)

:save frequency:
   how often the fitted weights are saved to an output file; (default: 1)

:test error:
   name of the file containing the misclassified test data; (default: error.dat)

:test input:
   name of the input file containing the test data;

:test output:
   name of the file containing the classification results of test data; (default: output.dat)

 
:train input:
   name of the input file containing the training data;

:weights input:
   name of input file containing the weights; (default: weights_old.dat)

:weights output:
   name of the output file containing the fitted weights; (default: weights.dat)






Training and test data
""""""""""""""""""""""""""

The format of the file containing the training or test data is given by the number of input and output nodes.  In general, there is a line for each data point which contains the values of the input functions and corresponding output values.  For classification, the output values are either 0 or 1. The number of columns is thus equal to the number of input + number of output nodes: e.g., if there are *x* input functions and *y* output values, then the first *x* columns will be the values of these input functions for a given configuration, and the next *y* columns the values of the corresponding output vector. 

In the example provided, the training and test data have 14 input functions and 5 output values as detailed in this `publication <https://doi.org/10.1103/PhysRevLett.123.245701>`_.


Weights file
"""""""""""""

The weights file contains the optimized weights and information concerning the architecture of the NN:

 - line 1: number of hidden layers
 - line 2: NN architecture; number of nodes in input, hidden, and output layers
 - line 3: translation and scaling for the input and output values
 - line 4: number of fitting iterations
 - line 5 and beyond: weights

Running the code 
^^^^^^^^^^^^^^^^^^^^

After the code has been compiled, it is executed by calling

.. code:: console 

   ./ml-classification.exe input.dat > screen.out

where ``input.dat`` is the name of the main input file. The code will run in parallel with the number of threads defined in the `input file <Main input file_>`_. 
In the command line above, the screen output is redirected to ``screen.out``.  

..
   In addition, the code writes two files, the ``weight.txt`` file and a file where the name is composed of the name of the input training data and the size of the hidden layers.  For example, if the input file is *Training.txt* and the NN has 2 hidden layers with 25 nodes each, the file is called ``Training_25_25_cost.txt``.

..
  :Training_25_25_cost.txt: 
   contains the cost function for each iteration; the meaning of the magnitude of the cost function is not entirely clear here, but it does need to converge to a constant value; this might not always be obvious to see

..
  :screen.out: 
   information about the architecture of the NN and running parameters, and again the convergence of the cost function

If the training is successful and the weights are converged, the file ``weights.dat`` (or whichever name was set for this file in the input) can be used to apply the NN for the classification of unknown environments.

Output files
"""""""""""""
In **training mode** (1 or -1), the code writes a weights file (default ``weights.dat``, see `Weights file`_), a cost function file (default ``cost.dat``), and prints some additional information to the screen. The cost function file contains the training loss and can be used to monitor the convergence of the fit.

In **testing mode** (0), the code writes an output and error file containing the predicted and true values of the test data.
The output of the NN is a continuous value :math:`x \in [0,1]` for each component (class) of the output layer.  For the test data, the true value is strictly either 0 or 1.  To compare the predicted values with the true values in the test file, the code checks if all components that should be 0 are :math:`< 0.5` AND if the component that should be 1 is :math:`> 0.5`.  If this is not the case, this data point is considered as wrongly predicted by the NN.

The L1 accuracy is defined here as 1-(wrong_points/total_points).

:screen.out: 
   standard output contains information about the architecture of the NN and running parameters; the L1 accuracy of the predicted classification with respect to the test set

:output.dat:
   for all test data points, the first *x* columns are the values of the *x* input functions, the next *y* columns are the predicted output values of the NN with *y* output nodes, and the last *y* columns are the actual values from the training set

:error.dat:
   same format as *output.txt*, but here only the points are recorded for which the NN failed to predict the correct class as defined above



Test examples
^^^^^^^^^^^^^^^^^^^^

To test the quality of a fitted NN for classification, the predicted values can be compared to known values for a given set of data points (e.g. different local structural environments).


In the folder ``testdata``, there are several input files for test data in various structural environments: ``test_BCC.txt``, ``test_FCC.txt``, ``test_HCP.txt``, ``test_A15.txt``, ``test_LIQ.txt``, and ``test_INT.txt`` within a body-centred cubic (BCC), face-centred cubic (FCC), hexagonal close-packed (HCP), A15, liquid (LIQ), and BCC/A15 interface (INT) structure (see also this `publication <https://doi.org/10.1103/PhysRevLett.123.245701>`_).

To run a test on any of these test data sets, provide a weight file, set the run option to 0, and set the name of the test input to any of the provided files.


.. warning:: 

   Make sure that the architecture (input, output, hidden layers) used with the given weights file is the same as used for training the neural network!



Compilation
^^^^^^^^^^^^^^^^^^^^

The NN fitting and testing code can be used either with `Intel MKL <https://software.intel.com/en-us/mkl>`_ or `OpenBLAS <https://www.openblas.net>`_ which needs to be installed on the system the code is supposed to run on.

(New York University cluster specific setup is explained :ref:`here <NYU>`)

Linux machines
"""""""""""""""
Using Intel MKL
++++++++++++++++++
Make sure `Intel MKL <https://software.intel.com/en-us/mkl>`_ is installed and environment variables (such as ``MKLROOT``) are set correctly.

Then the provided ``Makefile`` can be used and the code is compiled by

.. code:: console

   make

which will produce an executable ``ml-classification.exe``.  This executable will run in parallel with the number of threads given in the input file.


Using OpenBLAS
+++++++++++++++++++++
Make sure `OpenBLAS <https://www.openblas.net>`_ is installed and environment variables (such as ``BLAS_LIB``) are set correctly.

Then the provided ``Makefile`` can be used and the code is compiled by

.. code:: console

   make OBLAS=1

which will produce an executable ``ml-classification.exe``.  This executable will run in parallel with the number of threads given in the input file.


Mac
""""
Using OpenBLAS
+++++++++++++++++++++

If `OpenBLAS <https://www.openblas.net>`_ is not installed, this can, for example, be done using `MacPorts <https://www.macports.org>`_.  Usually, the corresponding ``INCLUDE`` and ``LDLIBS`` path should then be

.. code:: console
   
  /opt/local/include
  /opt/local/lib


If this is not the case, the ``Makefile`` needs to be adapted correspondingly.

The code can then be compiled using

.. code:: console
  
  make MAC_OBLAS=1

which will produce an executable ``Machine_Learning.exe``.  This executable will run in parallel with the number of threads given in the input file.


