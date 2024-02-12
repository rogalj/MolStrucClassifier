#include "header.h"
#include <string.h>
#include <ctype.h>
#include <algorithm>

#ifndef DAFED_H
#define DAFED_H

using namespace std;

//const int MAXNUMBEROFNEIGHBORS = 500;
const int MAXNUMBEROFNEIGHBORS = 600;
const int MAXNUMBEROFNEIGHBORS_1 = 10;
const int MAXATOMNUM = 100;
//const int MAXTRIO = MAXNUMBEROFNEIGHBORS*(MAXNUMBEROFNEIGHBORS-1)/2
const int nilvalue = 33333333;
const long int FACTORIALS[17] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000};
const int MAXBIN=10000;
const int MAX_NODES=51;
#define MY_NEIGHMASK 0x3FFFFFFF



//function declaration (usually in header file)
namespace DAFED{
class Cmolpoint {
	private:

	public: 
		Cmolpoint();
		virtual ~Cmolpoint();
		double neighbors[MAXNUMBEROFNEIGHBORS];
		int n_neighbors;
		int *atom_ID;
		int center;
		double COvec[3];	//used for proto type
		double NNvec[3];	//used for proto type
		//vector infomation here. same as parameter.symftype
		double **vecpoint;
		//vector_ID. should have same amount as parameter.symftype*2
		int *vector_ID;
		//symf result will be stored here
/**/		double G2vecNN;		//used for proto type
/**/		double G2vecCO;		//used for proto type
/**/		double G3vecNN;		//used for proto type
/**/		double G3vecCO;		//used for proto type
/**/		double point2;		//used for proto type	
/**/		double point3;		//used for proto type
		//new symf result here
		double symf;
		//distance from target to each neighbour
		double neighdist[MAXNUMBEROFNEIGHBORS];
		double diff[MAXNUMBEROFNEIGHBORS][3];
		double fcv[MAXNUMBEROFNEIGHBORS];
		double dfcvdx[MAXNUMBEROFNEIGHBORS][3];
};


class CAtom {
	private:
		//The positions of the particle
	public:
		CAtom();
		virtual ~CAtom();
		double pos[3];
		double COvec[3];
		double NNvec[3];
		double **vecpoint;
		int eletype;
		int moltype;
		int sa_mol_id;
		double G2vec;		//Delete sometime
		//derivative of point/point vector representation
		double d2point[3];
		double d3point[3];
		double dG2CO[3];
		double dG2NN[3];
		double G3vec;		//Delete sometime
		double dG3CO[3];
		double dG3NN[3];

		double fa[MAXNUMBEROFNEIGHBORS];
		double fb[MAXNUMBEROFNEIGHBORS];
		double costheta[MAXNUMBEROFNEIGHBORS][MAXNUMBEROFNEIGHBORS];
    		//JR: also save the distance in cartesian and spherical coordinates
		double n_diff[MAXATOMNUM][3];
    		double n_r[MAXNUMBEROFNEIGHBORS];
    		double n_phi[MAXNUMBEROFNEIGHBORS];
    		double n_theta[MAXNUMBEROFNEIGHBORS];

		int n_neighbors;

		//symmetry functions
		double sfgtype[9];
		double *sfg;
		int *ichunk;

		//value of NN output
		double *NNout;
		double *NNgrad;
		double *NNoutnorm;	//out_i/sum_i out_i
		double NNoutsum;	//sum_i out_i
		//derivative of global order parameter, dQ/dxi
		double dQ[3];
		//derivative of NN output wrt symmetry functions dq(i)/dsfg^k(i)
		double dqdsfg[13];

		//some more stuff for the local qlm
		double realmi[17], imgmi[17];
		double arealmi[17], aimgmi[17];
		double Lrealmi[3][17], Limgmi[3][17];	//save them for each L value to calculate derivatives
		double qL;				//atomic q value (q6,q7,q8...)
		double aqL;				//atomic average q value (aq6)
		double dQLM_ii[3][MAXNUMBEROFNEIGHBORS_1][3];  // [L][nn][x,y,z]
		double dQLM_ij[3][MAXNUMBEROFNEIGHBORS_1][3];  // [L][nn][x,y,z]
		
		//for the global Q6
		double d_QLM_pairs_real[13][3];
		double d_QLM_pairs_img[13][3];
		double d_Q6_pairs[3];  //derivative of Q6 over pairs
		double dreal_ij[MAXNUMBEROFNEIGHBORS_1][13][3];   //derivatives for certain rij [tj][mi][x,y,z]
		double dimg_ij[MAXNUMBEROFNEIGHBORS_1][13][3];   //derivatives for certain rij
		double d_Q6_ij[MAXNUMBEROFNEIGHBORS_1][3];    	// [nn][x,y,z]
		
		double d_surface; 	//DK: distance to surface.
		bool   second_mol;	//DK: Check if the molecule has same number of atom assigned in input file 
};

class QLMderiv{
	private:

	public:
		QLMderiv();
		virtual ~QLMderiv();

		double dreal_ii[MAXNUMBEROFNEIGHBORS_1][17][3];
		double dimg_ii[MAXNUMBEROFNEIGHBORS_1][17][3];
		double dreal_ij[MAXNUMBEROFNEIGHBORS_1][17][3];
		double dimg_ij[MAXNUMBEROFNEIGHBORS_1][17][3];
};


class CParameter{
	private:

	public:
		CParameter();
		virtual ~CParameter();
		//double neighbourdistance;
		double rmin0,rmax0;
		double rmin1,rmax1;
		double sf2_eta, sf2_Rs;  //for SF G2
		double sf3_kappa;    	//for SF G3
		//DK added
		double Rskappa;
		double eta;
		int flag;
		
		//////////////////////////	
		bool *pointflag;
		bool *g2g3flag; 
		int symftype;		//number of symmetry function type
		int *symftypeLst;	//type of symmetry function 

		/////////////////////////////////////
		//
		int npairs;
		int nsfg;
		int nsfgtot;
		int nstein;
		int nsfg2,nsfg3;  //how many symmetryfunctions of type 2 and 3
		int nsfg2CO,nsfg2NN,nsfg3CO,nsfg3NN,nsfg2point,nsfg3point;	//number of symmetry functions of g2 CO,NN, g3 CO,NN, point g2 and g3.
		int nnout;
		int COvectype[2];
		int NNvectype[2];
		//vectype store element type of each vector. 
		int *vectype;
		string INPUTfile;		

		double *RskappaLst;
		double *etaLst;

		//properties of extended variables
		int nex; 	//no. of extended variables
		int ex_seed;		//random number seed for extended variables
		double lwall; 		//lower wall for the sum of s1+s2 (here mainly A15+BCC)
		int lwall_cv;		//if wall potential acts on 0=exvar or 1=cv
		int lwall_n;		//parameters to determine wall potential: V(s)=k*((s-slimit)/eps)^n
		double lwall_k;
		double lwall_eps;
		
		int useNN;			// do we use a NN to get structure
		int choose_exvar;		// which exvar are we calculating: 0=QBCC/QA15, 1=f(QBCC,QA15), 2=Q6

		int use_restraint;		//additional restraint on sum BCC+A15, only works with f(QBCC,QA15)
		double restraint_value;		//where to place the restraint
		int restraint_n;		//restraining potential is V(s)=k*((s-slimit)/eps)^n
		double restraint_k;
		double restraint_eps;
		double restraint_pref;

		int useMetadyn;			// do we use a metadynamics bias
		int gauss_freq;			// how often a new gaussian is added
		double gauss_h, gauss_sigma;	// height and width of gaussians
		int bias_stride;		// how often the bias potential is printed
		int bias_read;			// if bias is read from file 0=no, 1=yes
		string bias_readfile;		// name of bias file

		int nmol;			//DK: number of molecules
		int center;			//DK: center of the mass carbon
		int natm;			//DK: number of atoms in the molecule
		bool triclinic;			//DK: added to check if the box is triclinic or not
		bool onlysymcal;
		bool NNrdf;			//DK: flag for surface distance calculation.
		int surface;			//DK: number of neighbor you use to determine surface
};

class CGgmt{
	// Second-order Generalized Gaussian Moments Thermostat
	private:

	public:
		CGgmt();
		virtual ~CGgmt();
		// Static parameters
		double Q1;	// thermostat mass1
		double Q2;	// thermostat mass2
		double kT0;	// reference temperature times the Boltzmann factor
		double dt[4];	// sub-time steps for internal RESPA and the Suzuki-Yoshida factorization
		int n_respa_ggmt;	//number of iterations for internal RESPA (default 1)

		//Dynamic variables
		double v1; 	// velocity1
		double v2;	// velocity2
		double eta1;	// position1
		double eta2;	// position2

};

class CHisto{
	// some variables to record a histogram
	private:

	public:
		CHisto();
		virtual ~CHisto();
		int nbin;			// no. of bins
		int nlat;			// no. of grid points (nlat = nbin+1)
		double min, max;		// minimum and maximum value of the histogram
		double binwidth;		// width of each bin
		int NconvDim;			// factor to convert from multi-dimensional array to 1D
						// e.g. in 3D we have a nlat_0 = 5, nlat_1 = 3, nlat_2 = 4 -> ntotlat = 5*3*4 = 60
						// NconvDim_0 = 60/5 = 12 , NconvDim_1 = 12/3 = 4 , NconvDim_3 = 4/4 = 1
						// (1,1,2) -> 1*12 + 1*4 + 2*1 = 18
						// 18 -> 18/12 = 1 and 18%12 = 6
						//        6/4  = 1 and 6%4   = 2
       						//        2/1  = 2		
};

//JR start: included for ramping s up/down
class CRamp{
	//some variables to ramp up/down the order parameter
	private:

	public:
		CRamp();
		virtual ~CRamp();
		double smin,smax;		//min and max value
		double ds;			//how much to increase/decrease s in each MD step
		int steps;			//how many steps between smin and smax steps = (smax-smin)/ds
		int steps0;			//initial no. of steps from current s to smax
		int count;			//counter to check how many steps have been done
		double switchme;		//+1/-1 depending on where to move s
};
//JR end: included for ramping s up/down

class EXvariables{
	private:
	
	public:
		EXvariables();
		virtual ~EXvariables();
		double kappa, tau, gamma, temp;
		double Q,Qlocal;		// value of cv, Qlocal is value owend by certain proc
		double QU1,Qliq;		//DK: This is for testing Q
		double **dQ;			// derivative wrt atom i
		double **dQU1;			// DK: This is for test. derivative wrt atom i
		double **dQUliq;		// DK: This is for test. derivative wrt atom i
		double *qeach;
		double **dqeach;

		double x,v,f;			// position, velocity, total force
		double fcoup, fbias;		// harmonic force from coupling and force from additional bias
		double fwall;			// force from lower wall of sum
		double fLwall,fUwall;		// force from a lower/upper bound
		double fatom;			// actual force that is added to atomic forces!
		double frestrain;		// additional force due to restraining potential
		double vrestrain;		// value of the restraining potential
		double m;			// mass
		double dvirial[3][3];		// derivative part of the virial \sum dq/dx*x [x][fx]
		double c1_lgvn, c2_lgvn;	// coefficients for langevin thermostat
		CGgmt ggmt;			// variables for ggmt thermostat

		int n_respa;			// Number of RESPA integration steps between two MD steps
		double dt_respa;		// sub-time step for DAEFE RESPA

		double dafed_work;		// work of meta_variable on physical system
		double dWold;			// Old work element for trapeze integration
		double thermo_work;		// Thermostat work

		double lwall,uwall;		// values of lower and upper wall, slimit
		int lwall_n, uwall_n;		//parameters to determine wall potential: V(s)=k*((s-slimit)/eps)^n
		double lwall_eps, lwall_k;
		double uwall_eps, uwall_k;
		double uwall_pref, lwall_pref;	//prefactor = -k*n/eps^n
		int lwall_cv, uwall_cv;		//if constraint is applied to 0=exvar or 1=cv
		
		CHisto histo;			//parameters for the bias histogram
		
		//JR start: included for ramping s up/down
		CRamp ramp;			//parameters to ramp up/down s
		//JR end: included for ramping s up/down
};



class NNvariables{ 
	private: 

	public:
		NNvariables();
		virtual ~NNvariables();

		int layers; 	//amount of hidden layers
		int input;  	//amount of inputs
		int output;	//amount of outputs
		int *arch;	//amount of nodes in layer n

		int activation;  // which activation function to use; 0=sigmoid, 1=softmax; also changes loss function accordingly 
		
		double *w;  	//weights in the synapses
		double *nodes;	//values in each node after networking once (a terms), h(a)=z terms
		double *Gnodes;	//y values in each node after networking once (a terms), h(a)=z terms
		double *translation;  	//translation for feature scaling
		double *stretch;  	//stretch for feature scaling

}; 

//for neighbour lists from lammps
class CNeighvariables{  
	private:

	public:
		CNeighvariables();
		virtual ~CNeighvariables();

		int inum;				// # of neighbour lists in lammps
		int *ilist, *numneigh, **firstneigh;	// other stuff of neighbour lists
		int nall;				// no of owned+ghost atoms
		int ntot;				// total number of atoms in cell
		int lmpmol;			//DK: added for lmpmol
		
		int me;		//only for now the proc no.
};


struct ComplexNumber{
	double real,img;
};

//this is for the histogram of the bias potential in CV space
struct Bhist{
	double *x;		// value of the CV, x is array of length N_CV
	double pot;		// value of the bias potential
};

//this is for the interpolation of the bias potential in between grid points, 2^d vertices in d-dimensional space
struct Nblat{
	double sign;		// +1/-1 depending on the vertex
	int *diff;		// vector that contain 0/1 for each CV depending on the vertex
};




	//reading dafed input file
	void read_dafed_input(const string &, CParameter &, EXvariables *, string &, int &, string &, int &);

	//random number generator
	double ran3(int=1);
	double gauss();

	//test symmetry function genreralized code
	void symf_multi_vec_sa(Cmolpoint *, CParameter &, int, int); //typeNum from symfType. sfgNum from i in nsfg loop
	void symf_multi_vec_sa_short_test(CAtom *, Cmolpoint *, CParameter &); //typeNum from symfType. sfgNum from i in nsfg loop
	void symf_multi_vec_lmp_short_test(CAtom *, Cmolpoint *, CParameter &,CNeighvariables &); //typeNum from symfType. sfgNum from i in nsfg loop
	void symf_multi_vec_lmp(CAtom *, Cmolpoint *, CParameter &, CNeighvariables &, int, int);
//call the each symmetry function in here
	void get_vec_symmetry_NN_3_lmpneighbour_new(CAtom *, Cmolpoint *,CNeighvariables &, double *, CParameter &, double *, double *,  NNvariables &);
	void get_vec_symmetry_NN_vec_sa(CAtom *, Cmolpoint *,int &, CParameter &,  NNvariables &);
	void get_vec_symmetry_NN_vec_sa_short_test(CAtom *, Cmolpoint *,int &, CParameter &,  NNvariables &, ofstream &);
	void get_vec_symmetry_NN_vec_lmp_short_test(CAtom *, Cmolpoint *, CNeighvariables &lmpneigh, CParameter &,  NNvariables &);
	void get_vec_symmetry_NN_3_lmpneighbour_vec_lmp(CAtom *, Cmolpoint *, CNeighvariables &, CParameter &, NNvariables &);
	void get_vec_symmetry_output(CAtom *, Cmolpoint *,int &, CParameter &, ofstream &);
	void get_NNrdf(CAtom *, Cmolpoint *,int &, CParameter &, double *);
	void get_print_NNrdf(CAtom *, Cmolpoint *,int &, CParameter &, NNvariables &, ofstream &);
	void print_traj(CAtom *, Cmolpoint *, double *, int &, double **, CParameter &, NNvariables &, double *, int &, int &,int *,ofstream&, ofstream&, ofstream&, ofstream&);


//find the neighbor information in here	
	void get_allNeighbourDistances(CAtom *, Cmolpoint *, int &, double *, CParameter &);


	void get_allNeighbourDistances_lmp_new(CAtom *, Cmolpoint *, CNeighvariables &, int &, double *, CParameter &);
	void get_allNeighbourDistances_gen_vec_lmp(CAtom *, Cmolpoint *, CNeighvariables &, int &, double *, CParameter &);
	void get_allNeighbourDistances_sa(CAtom *, Cmolpoint *, int &, double *, CParameter &);
	void get_allNeighbourDistances_gen_vec_sa(CAtom *, Cmolpoint *, int &, double *, CParameter &);

	double get_absDistance(int ,int ,double & ,double &,double &, CAtom *, double *, bool);
	void get_molvec_pbc(int,int,double *, CAtom *, double *, bool);
	//DK added
	//DK: This calculates each symmetry functions 
	void symmetryfunc_g2_fc0_vecCO_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);	
	void symmetryfunc_g3_fc0_vecCO_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);
	void symmetryfunc_g3_fc0_vecNN_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);
	void symmetryfunc_g2_fc0_vecNN_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);
	
	void symmetryfunc_g2_fc0_vecCO_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symmetryfunc_g2_fc0_vecNN_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symmetryfunc_g3_fc0_vecCO_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symmetryfunc_g3_fc0_vecNN_sa(CAtom *, Cmolpoint *, int &, CParameter &);

	//DK: added for lmp neighbor
	void symfg2_lmp(CAtom *, Cmolpoint *, int &, CParameter &,CNeighvariables &);
	void symfg3_lmp(CAtom *, Cmolpoint *, int &, CParameter &,CNeighvariables &);
	
	void symfg2_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symfg3_sa(CAtom *, Cmolpoint *, int &, CParameter &);

	//calculate the derivative of each symmetry function	
	void get_total_derivatives_molvec_lmp(CAtom *, Cmolpoint *, int &,CParameter &, int *, int *, EXvariables *, double *, double *,CNeighvariables &);
	void get_total_derivatives_molvec_sa(CAtom *, Cmolpoint *, int &,CParameter &, EXvariables *, double *, double *);
	void get_total_derivatives_molvec_vec_sa(CAtom *, Cmolpoint *, int &,CParameter &, EXvariables *);
	void get_total_derivatives_molvec_vec_lmp(CAtom *, Cmolpoint *, int &,CParameter &, EXvariables *,CNeighvariables &);

	
	void get_Q_vec_lmpneigh_vec_new(CAtom *, Cmolpoint *, CNeighvariables &, EXvariables *, CParameter &, double *);

	//NN functions
	void read_weight( NNvariables &);
	void get_NN(CAtom *,int, int &, NNvariables &);  
	void get_NN_norm(CAtom *,int &, NNvariables &);  
	void get_NN_test(double **,int &, double **, double **, NNvariables &);  
	void get_NN_vec(Cmolpoint *,int &, NNvariables &);
	
	
	void NN_network(double *, double * , double *, NNvariables &);
	void NN_hidden( NNvariables &, int ); 
	void NN_step(NNvariables &, int , double *);
	void NN_stepGradient(NNvariables &, int , double *);
	double h_class(double );
	double hprime_class(double );
	double softmax(double, double);
	double softmaxprime(double, double, int, int, double);




	//only for  testing
	//void test_aq6(double **,int &, double *,CParameter &);
	void test_vec(double **,int &, int *, int *, double *,CParameter &,string &,double *,int &, int *);
	void read_input(const string &, CParameter &);


}

#endif // DAFED_H
