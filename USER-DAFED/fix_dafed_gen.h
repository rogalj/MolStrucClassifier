
#ifdef FIX_CLASS

FixStyle(dafedgen,FixDafedgen)

#else

#ifndef LMP_FIX_DAFEDGEN_H
#define LMP_FIX_DAFEDGEN_H

#include "fix.h"
#include "compute.h"
#include <fstream>
#include <sstream>

//#include <stdlib.h>
//#include <cstdlib>
//#include <stdio.h>

// the dafed header that defines the class//
#include "../Dafed.h"

using namespace std;

namespace LAMMPS_NS {

class FixDafedgen : public Fix {
 	public:
  	FixDafedgen(class LAMMPS *, int, char **);
  	~FixDafedgen();
  	int setmask();
  	void init();
	void init_list(int, class NeighList *);
	void setup(int);
  	void post_force(int);
  	void post_force_respa(int, int, int);
	void end_of_step();
	
	void integrate_ggmt(DAFED::EXvariables &);

	void init_colvar_file(string &);			//open colvar file and write header lines
	void init_nnout_file();					//open nnout file and write header line
	void init_dafedene_file();				//open dafed_ene file and write header lines
	void init_bias_files();					//open files to record forces and bias potential on extended variable in UFED
	void init_derivatives_files(string &);				//open files and write headers for derivatives, one for each proc

	void print_colvar(int, double);			//print colvar/exvar values to file	
	void print_nnout(int, double);			//print NN output to file
	void print_dafed(DAFED::EXvariables &);		//print dafed energy and work transferred to the system
	void print_biasforce(DAFED::EXvariables &);	//print force due to different biases (gaussians and upper/lower walls)
	void print_biaspot();				//print the gaussian bias potential
	
	void get_exvar_values();
	void get_pathcv_values();
	void get_pathcv_values_Xpoints();

	void get_pathcv_values_test();			//DK: Just for test
	void get_pathcv_values_gen();
	void get_Q6_values();
	void get_numerical_derivatives();

	void setup_ufedbias();			// setup histograms etc to record additional metadynamics bias
	void update_ufedbias();			// update Gaussian bias
	void read_ufedbias();			// read bias potential from file
	void get_bias_force();			// get force from the bias potential
	void get_bias_fixedforce();		// set a fixed value of the bias force
	void get_ind_rev(int *, int); 		// convert 1D array index into index for each exvar
	int get_ind(int *);			// convert n-dim index for all exvar into 1D index

	//stuff from fcc/orient
	struct Nbr {              	// neighbor info for each owned and ghost atom
    		tagint id[30];          // IDs of each neighbor
                            		// if center atom is owned, these are local IDs
                            		// if center atom is ghost, these are global IDs
  	};
	int pack_forward_comm(int, int *, double *, int, int *);
  	void unpack_forward_comm(int, int, double *);
	int pack_reverse_comm(int, int, double *);
	void unpack_reverse_comm(int, int *, double *);
	
	double memory_usage();

	void grow_arrays(int);			//for growing atom based arrays
	void copy_arrays(int, int, int);	//to copy from atom i to j; for sorting?
	int pack_exchange(int, double *);	//store per atoms values when atoms migrate...
	int unpack_exchange(int, double *);


 	private:
	// output file
	ofstream of;
	ofstream of_colvar;
	ofstream of_deriv;
	ofstream of_nnout;
	ofstream of_ene;
	ofstream of_biasforce;
	ofstream of_biaspot;
	//begin debug
	//ofstream of_debug_virial;
	//end debug
	
	// this is something to enable respa
  	int nlevels_respa;
	
	//how often to print someting
	int nprint_cv,count_cv,count2_cv;
	int nprint_deriv,count_deriv;
	int nprint_bias,count_printbias;

	// Compute for the energy
	class Compute *c_pe; 

	double **dQ6_pairs; //derivatives of Q6_pairs
	double Q6_pairs;

	//Parameters for cutoff/symmetry/...
	DAFED::CParameter parameter;
	//variables for the NN;
	DAFED::NNvariables NNvar;
	//extended variables class 	
	DAFED::EXvariables exvar[2];
	DAFED::CNeighvariables lmpneigh;
	//parameters for the symmetry functions
	//when using 13 symmetry functions
	//double sf2Rs[13]={2.3, 2.4, 2.5, 2.8, 3.1, 3.5, 4.0, 4.4, 4.4, 4.5, 5.1, 5.2, 5.4};
	//double sf2eta[13]={20.0, 0.5, 0.5, 0.5, 0.5, 0.5, 5.0, 0.5, 20.0, 0.5, 5.0, 2.5, 5.0};
	//when using 8 symmetry functions
	//double sf2Rs[8]={2.3, 2.4, 4.0, 4.4, 4.4, 5.1, 5.2, 5.4};
	//double sf2eta[8]={20.0, 0.5, 5.0, 0.5, 20.0, 5.0, 2.5, 5.0};

	//for 11 symmetry functions, 8 sfg2 and 3 sfg3
	//CO G2G3 0,1
        //NN G2G3 2,3
        //Sy G2G3 4,5              
	double Rskappa[24]={6.16,6.28,6.76,6.88,0.36,0.08,0.36,0.28,-0.64,-0.36,0.88,1.0,2.5,4.54,4.9,6.22,2.5,3.58,4.78,8.26,2.50,8.12,8.24,8.36};	
	double eta[24]={2.44,2.68,1.0,1.0,1.0,1.0,1.12,6.76,3.28,3.28,3.28,3.28,1,1,1,1,1,1,1,1,1,1,1,1};
	
	int totalsymf= -1;		//number of symf you use for this calculation. Initialized as -1
	double *RskappaLst;	
	double *etaLst;
	int *sfgnum;
	//int Calcflag[24] = {4,4,4,4,0,0,0,0,2,2,2,2,5,5,5,5,1,1,1,1,3,3,3,3};

	//double sf3kappa[3]={3.5, 4.5, 7.0};


	//stuff for neighbor lists
	class NeighList *list;
	DAFED::CAtom *molecules;		//data structure used within the library
	DAFED::Cmolpoint *molec;
	//void get_atominfo(int, CNeighvariables &, CAtom *);
	void get_atominfo(int);
	double Qall[6];
	double Qlocal[6];

	

	//stuff from orient_fcc
	Nbr *nbr;
	int nmax;		//max number of owned + ghost atoms on this proc
	//DK:added for communication
	int commflag;
	int me;			//which proc is this
	int nprocs;		//total number of procs

	//stuff for recording and interpolating bias potential
	int Nvertex;		// no. of vertices for d-dimensional interpolation, Nvertex = 2^d
	int Nlattot;		// total number of lattice points
	DAFED::Bhist *bhist;		// 1D array that contains all the lattice points; length is \prod n_i , where n_i is grid point in direction i
				// each lattice point has a value 'pot' and a vector 'x'
	DAFED::Nblat *nblat;		// this is the neighbour lattice for the interpolation between grid points
	int count_updatebias;

	//for a fixed bias force (no bias potential required)
	double biasforce_sign;

	// array to store CV derivatives for each atom
	double **array;
	int nvalues;

	//DK:added for printing out NNout
	double **NNout;
};

};

#endif
#endif
