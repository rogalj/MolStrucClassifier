//#include "header.h"
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
//		int atom_ID[4];
/**///		int atom_ID[8];
//		int atom_ID[15];
		int *atom_ID;
		int center;
		double COvec[3];
		double NNvec[3];
		double G2vecNN;
		double G2vecCO;
		double G3vecNN;
		double G3vecCO;
		double point2;
		double point3;
		double neighdist[MAXNUMBEROFNEIGHBORS];
		double diff[MAXNUMBEROFNEIGHBORS][3];
/**/		double sfg[24];
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
		int eletype;
		int moltype;
		double G2vec;		//Delete sometime
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
		//double ac, Rc, ec;  	//for fc
		//double eta, mu; 	//for fa
		//double nu, al, ar;	//for fb
		double sf2_eta, sf2_Rs;  //for SF G2
		double sf3_kappa;    	//for SF G3
		//DK added
		double Rskappa;
		double eta;
		int flag;

		//double sf4_eta, sf4_lambda, sf4_zeta;	//for SF G4
		//double sf5_eta, sf5_lambda, sf5_zeta;	//for SF G5
		//double sf6_lambda, sf6_zeta;		//for SF G6
		//double sf7_eta, sf7_alpha;		//for SF G7
		//double sf8_eta, sf8_alpha;		//for SF G8
		int npairs;
		int nsfg;
		int nsfgtot;
		int nstein;
		int nsfg2,nsfg3;  //how many symmetryfunctions of type 2 and 3
		int nsfg2CO,nsfg2NN,nsfg3CO,nsfg3NN,nsfg2point,nsfg3point;	//number of symmetry functions of g2 CO,NN, g3 CO,NN, point g2 and g3.
		int nnout;
		int COvectype[2];
		int NNvectype[2];
		double *RskappaLst;
		double *etaLst;

		//properties of extended variables
		int nex; 	//no. of extended variables
		//double ex_kappa[5],ex_tau[5],ex_gamma[5],ex_temp[5];	//kappa, tau, gamma, and T of extended variable
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



	//Functions for symmetry function
//	void get_QBCC(double **, int &, double *, double **, double **, CParameter &, double *, double *, NNvariables &);
//	void get_QBCC_short(double **, int &, double *, double **, double **, CParameter &, double *, double *, NNvariables &);
//	void get_Q_short(double **, int &, double *, EXvariables *, CParameter &, double *, double *, double *, NNvariables &, double *);
//	void get_Q_short_3(double **, int &, double *, EXvariables *, CParameter &, double *, double *, double *, NNvariables &, double *);
//	void get_Q_3_lmpneighbour(CAtom *, CNeighvariables &, EXvariables *, CParameter &, double *, double *, double *, NNvariables &, double *);
	
	//this is also in steinhardt.cpp
//	void calculate_globalQ6_lmpneighbour_part01(CAtom *, int &, int &, double *, double *);
//	void calculate_globalQ6_lmpneighbour_part02(CAtom *, int &, int &, double *, double *, EXvariables *);
//	void calculate_globalQ6_virial_lmpneighbour_part01(CAtom *, int &, int &, double *, double *);
//	void calculate_globalQ6_virial_lmpneighbour_part02(CAtom *, int &, int &, double *, double *, EXvariables *);
	
//	void get_Q_short_stein(double **, int &, double *, EXvariables *, CParameter &, NNvariables &, double *);
//	void get_Qnorm_short(double **, int &, double *, EXvariables *, CParameter &, double *, double *, double *, NNvariables &, double *);
	
//	void get_total_derivatives(CAtom *, int &, CParameter &, double *, double *);
//	void get_total_derivatives_short(CAtom *, int &, CParameter &, double *, double *, EXvariables *);
//	void get_total_derivatives_short_2(CAtom *, int &, CParameter &, double *, double *, double *, EXvariables *);
//	void get_total_derivatives_short_3(CAtom *, int &, CParameter &, double *, double *, double *, EXvariables *);
//	void get_total_derivatives_3_lmpneigh(CAtom *, int &, int &, CParameter &, double *, double *, double *, EXvariables *);
//	void get_total_derivatives_short_stein(CAtom *, int &, CParameter &, EXvariables *);
	//void get_Q_derivatives(CAtom *, int &, CParameter &, double *, double *, double *, EXvariables *);
	
//	void get_total_derivatives_short_norm(CAtom *, int &, CParameter &, double *, double *, EXvariables *);
//	void get_total_derivatives_short_norm_2(CAtom *, int &, CParameter &, double *, double *, double *, EXvariables *);
	
//	void get_symmetry(CAtom *, int &, double *, double **, CParameter &, double *, double *);
//	void get_symmetry_short(CAtom *, int &, double *, CParameter &, double *, double *);
//	void get_symmetry_short_2(CAtom *, int &, double *, CParameter &, double *, double *, double *);
//	void get_symmetry_short_3(CAtom *, int &, double *, CParameter &, double *, double *, double *);
//	void get_symmetry_NN_3_lmpneighbour(CAtom *, CNeighvariables &, double *, CParameter &, double *, double *, double *, NNvariables &);
//	void get_vec_symmetry_NN_3_lmpneighbour(CAtom *, Cmolpoint *,CNeighvariables &, double *, CParameter &, double *, double *,  NNvariables &, EXvariables *);
	void get_vec_symmetry_NN_3_lmpneighbour_new(CAtom *, Cmolpoint *,CNeighvariables &, double *, CParameter &, double *, double *,  NNvariables &);
	void get_vec_symmetry_NN_sa(CAtom *, Cmolpoint *,int &, CParameter &, double *, double *,  NNvariables &);
//	void symfg2(CAtom *, Cmolpoint *, int &, CParameter &);
//	void symfg3(CAtom *, Cmolpoint *, int &, CParameter &);

//	void get_q6_lmpneighbour(CAtom *, CNeighvariables &, double *, CParameter &);
//	void get_aq6_lmpneighbour(CAtom *, CNeighvariables &, CParameter &, double &, int &);

//	void get_symmetry_short_stein(CAtom *, int &, double *, CParameter &);
	
	void get_allNeighbourDistances(CAtom *, Cmolpoint *, int &, double *, CParameter &);
//	void get_allNeighbourDistances_cutoff(CAtom *, int &, double *, CParameter &);
//	void get_allNeighbourDistances_2Xcutoff(CAtom *, int &, double *, CParameter &);

//	void get_info();

//	void get_allNeighbourDistances_cutoff_short(CAtom *, int &, double *, CParameter &);
	//void get_allNeighbourDistances_cutoff_lmpneighbour(CAtom *, CNeighvariables &, double *, CParameter &);
//	void get_allNeighbourDistances_lmp(CAtom *, Cmolpoint *, CNeighvariables &, int &, double *, CParameter &);
	void get_allNeighbourDistances_lmp_new(CAtom *, Cmolpoint *, CNeighvariables &, int &, double *, CParameter &);
	void get_allNeighbourDistances_sa(CAtom *, Cmolpoint *, int &, double *, CParameter &);
	//void TESTfunc(CAtom *, Cmolpoint *, CNeighvariables &, int &, double *, CParameter &);
//	void TESTfunc(CAtom *, Cmolpoint *, CNeighvariables &, double *, CParameter &, double *, double *,  NNvariables &);
//	void get_allNeighbourAngles(CAtom *, int &);
	double get_absDistance(int ,int ,double & ,double &,double &, CAtom *, double *, bool);
	void get_molvec_pbc(int,int,double *, CAtom *, double *);
	//void get_cutofffunctions(CAtom *, int &, CParameter &);
	//void symmetryfunc_g1(CAtom *, int &);
	//void symmetryfunc_g2(CAtom *, int &, CParameter &);
//	void symmetryfunc_g2_fc0(CAtom *, int &, CParameter &);
//	void symmetryfunc_g2_fc0_short(CAtom *, int &, CParameter &);
//	void symmetryfunc_g3_fc0_short(CAtom *, int &, CParameter &);
//	void symmetryfunc_g3_fc0(CAtom *, int &, CParameter &);
	//DK added
//	void symmetryfunc_g2_fc0_vecCO(CAtom *, int &, CParameter &, int *);
//	void symmetryfunc_g2_fc0_vecCO_newTEST(CAtom *, Cmolpoint *, int &, CParameter &, int *, int &,EXvariables *);	
//	void symmetryfunc_g3_fc0_vecCOTEST(CAtom *, Cmolpoint *, int &, CParameter &, int *, int &, EXvariables *);
//	void symmetryfunc_g3_fc0_vecNNTEST(CAtom *, Cmolpoint *, int &, CParameter &, int *, int &, EXvariables *);
//	void symmetryfunc_g2_fc0_vecNNTEST(CAtom *, Cmolpoint *, int &, CParameter &, int *, int &, EXvariables *);
	//DK: added for lmp neighbor
	void symmetryfunc_g2_fc0_vecCO_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);	
	void symmetryfunc_g3_fc0_vecCO_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);
	void symmetryfunc_g3_fc0_vecNN_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);
	void symmetryfunc_g2_fc0_vecNN_lmp(CAtom *, Cmolpoint *, int &, CParameter &, int *,CNeighvariables &);
	
	void symmetryfunc_g2_fc0_vecCO_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symmetryfunc_g2_fc0_vecNN_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symmetryfunc_g3_fc0_vecCO_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symmetryfunc_g3_fc0_vecNN_sa(CAtom *, Cmolpoint *, int &, CParameter &);

	//void symmetryfunc_g2_fc0_vec_genTEST(CAtom *, Cmolpoint *,int &, CParameter &, int *,int &, EXvariables *);
	//void symmetryfunc_g3_fc0_vec_genTEST(CAtom *, Cmolpoint *,int &, CParameter &, int *,int &, EXvariables *);	
//	void symfg2TEST(CAtom *, Cmolpoint *, int &, CParameter &, int &, EXvariables *);
//	void symfg3TEST(CAtom *, Cmolpoint *, int &, CParameter &, int &, EXvariables *);
	//DK: added for lmp neighbor
	void symfg2_lmp(CAtom *, Cmolpoint *, int &, CParameter &,CNeighvariables &);
	void symfg3_lmp(CAtom *, Cmolpoint *, int &, CParameter &,CNeighvariables &);
	
	void symfg2_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	void symfg3_sa(CAtom *, Cmolpoint *, int &, CParameter &);
	
//	void symmetryfunc_g2_fc0_vecCO_new(CAtom *, Cmolpoint *, int &, CParameter &, int *);	
//	void symmetryfunc_g3_fc0_vecCO(CAtom *, Cmolpoint *, int &, CParameter &, int *);
//	void symmetryfunc_g3_fc0_vecNN(CAtom *, Cmolpoint *, int &, CParameter &, int *);
//	void symmetryfunc_g2_fc0_vecNN(CAtom *, Cmolpoint *, int &, CParameter &, int *);
//	void symmetryfunc_g2_fc0_vec_gen(CAtom *, Cmolpoint *,int &, CParameter &, int *);
//	void symmetryfunc_g3_fc0_vec_gen(CAtom *, Cmolpoint *,int &, CParameter &, int *);	
//	void get_total_derivatives_molvecg2(CAtom *, int &, CParameter &, int *);
//	void get_total_derivatives_molvecg2_CO(CAtom *, Cmolpoint *, int &, CParameter &, int *);
//	void get_total_derivatives_molvecg2_NN(CAtom *, Cmolpoint *, int &, CParameter &, int *);
//	void get_total_derivatives_molvecg3_CO(CAtom *, Cmolpoint *, int &,CParameter &, int *);
//	void get_total_derivatives_molvecg3_NN(CAtom *, Cmolpoint *, int &,CParameter &, int *);
//	void get_total_derivatives_molvecg3(CAtom *, int &, CParameter &, int *);
//	void get_total_derivatives_pointg2(CAtom *, Cmolpoint *, int &,CParameter &);
//	void get_total_derivatives_pointg3(CAtom *, Cmolpoint *, int &,CParameter &);
	
//	void get_vec_allNeighbourDistances_cutoff_lmpneighbour(CAtom *, Cmolpoint *, CNeighvariables &, double *, CParameter &);
	
//	void get_total_derivatives_molvec(CAtom *, Cmolpoint *, int &,CParameter &, int *, int *, EXvariables *, double *, double *);
	void get_total_derivatives_molvec_lmp(CAtom *, Cmolpoint *, int &,CParameter &, int *, int *, EXvariables *, double *, double *,CNeighvariables &);
	void get_total_derivatives_molvec_sa(CAtom *, Cmolpoint *, int &,CParameter &, EXvariables *, double *, double *);
	void get_Q_vec_lmpneigh(CAtom *, Cmolpoint *, CNeighvariables &, EXvariables *, CParameter &, double *, double *, NNvariables &, double *);
	//void symmetryfunc_g4(CAtom *, int &, CParameter &);
	//void symmetryfunc_g5(CAtom *, int &, CParameter &);
	//void symmetryfunc_g6(CAtom *, int &, CParameter &);
	//void symmetryfunc_g7(CAtom *, int &, CParameter &);
	//void symmetryfunc_g8(CAtom *, int &, CParameter &);


	//Steinhardt parameter functions
//	void get_qL(CAtom *, int &, int &); 	//this is over pairs
//	void get_qL_full_value(CAtom *, int &, int &);	//this is over all neighbours for each atom, only value
//	void get_qL_full_derivative(CAtom *, int &, int &);	//this is over all neighbours for each atom, only derivative
	
//	void get_aqL_full_value(CAtom *, int &, int &);	//this is over all neighbours for each atom, only value
	
//	void QLM(int ,int ,double, double *, double &, double &, double *, double *);
//	void QLM_value(int ,int ,double, double *, double &, double &);
	
//	void YLM_nocomplex(int, int, double, double *, double &, double &, double *, double *);	
//	void YLM_nocomplex_value(int, int, double, double *, double &, double &);	
	
//	double PLM_new(int, int, double, double &);
//	double PLM_new_value(int, int, double);



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




	//only for  testing
	//void test_aq6(double **,int &, double *,CParameter &);
	void test_vec(double **,int &, int *, int *, double *,CParameter &);
	void read_input(const string &, CParameter &);


}

#endif // DAFED_H
