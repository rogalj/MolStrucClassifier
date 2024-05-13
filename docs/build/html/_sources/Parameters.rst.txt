Parameters
==========



Parameters you have to choose are in ``dAFED.in``, ``fix_dafed.cpp``,``compute_NNout.cpp``, and run file. 

In fix_dafed.cpp, check parameters around line 494

.. code-block:: console

   parameter.rmin0 = 9.8;
   parameter.rmax0 = 10;
   parameter.npairs = 0;
   parameter.nsfg=24;
   parameter.nsfg2CO=4;    
   parameter.nsfg2NN=4;    
   parameter.nsfg3CO=4;    
   parameter.nsfg3NN=4;    
   parameter.nsfg2point=4; 
   parameter.nsfg3point=4; 
   parameter.nmol = natoms/parameter.natm; 
   parameter.center = 3;
   parameter.nnout = 6;
   parameter.COvectype[0] = 3;
   parameter.COvectype[1] = 2;
   parameter.NNvectype[0] = 4;
   parameter.NNvectype[1] = 5;

In fix_dafed.h, check parameters around line 135

.. code-block:: console

   double Rskappa[24]={6.16,6.28,6.76,6.88,0.36,0.08,0.36,0.28,-0.64,-0.36,0.88,1.0,2.5,4.54,4.9,6.22,2.5,3.58,4.78,8.26,2.50,8.12,8.24,8.36};
   double eta[24]={2.44,2.68,1.0,1.0,1.0,1.0,1.12,6.76,3.28,3.28,3.28,3.28,1,1,1,1,1,1,1,1,1,1,1,1};


In compute_NNout.cpp, check parameters around line 318

.. code-block:: console

   parameter.rmin0 = cutoff_user-0.2;
   parameter.rmax0 = cutoff_user;
   parameter.nsfg = 24;
   parameter.center = 3;
   parameter.natm = 8;
   parameter.nsfg=24;
   parameter.nsfg2CO=4;    
   parameter.nsfg2NN=4;    
   parameter.nsfg3CO=4;    
   parameter.nsfg3NN=4;    
   parameter.nsfg2point=4; 
   parameter.nsfg3point=4; 
   parameter.nmol = natoms/parameter.natm; //no. of total mol in the system
   parameter.center = 3;
   parameter.nex = 2;
 


In compute_NNout.h, check parameters around line 77, 78. 

.. code-block:: console

   double Rskappa[24]={6.16,6.28,6.76,6.88,0.36,0.08,0.36,0.28,-0.64,-0.36,0.88,1.0,2.5,4.54,4.9,6.22,2.5,3.58,4.78,8.26,2.50,8.12,8.24,8.36};
   double eta[24]={2.44,2.68,1.0,1.0,1.0,1.0,1.12,6.76,3.28,3.28,3.28,3.28,1,1,1,1,1,1,1,1,1,1,1,1};


