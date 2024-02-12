#include "math.h"
//#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "fix_dafed_gen.h"
#include "universe.h"
#include "compute.h"
#include "modify.h"
#include "pair.h"
#include <iostream>
#include <iomanip>

#include <chrono>
#include <ctime>

#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

#include "comm.h"
//from fcc/orient
#include "memory.h"
#include <mpi.h>



using namespace LAMMPS_NS;
//using namespace PLMD;
using namespace FixConst;
using namespace std;

#define INVOKED_SCALAR 1

FixDafedgen::FixDafedgen(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  NNout(nullptr),
  me(0),
  nprocs(0),
  of(NULL),	
  of_colvar(NULL),
  of_deriv(NULL),
  of_nnout(NULL),
  of_ene(NULL),
  of_biasforce(NULL),
  of_biaspot(NULL),
  nprint_cv(0),
  nprint_deriv(0),
  nprint_bias(0),
  count_cv(0),
  count2_cv(0),
  count_deriv(0),
  count_printbias(0),
  //x_s(NULL),
  //v_s(0.0),
  //f_s(0.0),
  //m_s(0.0),
  //T_s(0.0),
  //kappa_s(0.0),
  //c1_lgvn(0.0),
  //c2_lgvn(0.0),
  parameter(),
  exvar(),
  list(NULL),
  nbr(NULL),	
  lmpneigh(),
  Qall(),
  Qlocal(),
  Nvertex(0),
  Nlattot(0),
  bhist(NULL),
  nblat(NULL),
  count_updatebias(0),
  biasforce_sign(1.0),
  //begin debug
  //of_debug_virial(NULL),
  //end debug
  molecules(),
  molec(),
  array(nullptr),
  nvalues(0)
{
	
	//This is necessary to properly communicate stuff:
  	if (!atom->tag_enable) error->all(FLERR,"fix dafed requires atom tags");

	//check if this is run on more than one core
	MPI_Comm_rank(world,&me);
  	MPI_Comm_size(world,&nprocs);
	if(nprocs > 1){
		cout<<"Hello, I'm running dAFED in parallel with nprocs : "<<nprocs<<endl;
		cout<<"This is proc no "<<me
			<<" nlocal = "<<atom->nlocal
			<<endl;
		//error->all(FLERR,"fix dafed only works on SINGLE CORES right now!");
	}

	//from fcc/orient, initialize some stuff for mpi communication
	peratom_flag = 1;		//calculate/communciate per atom quantities
	nmax=0;				//for array allocation later
	//I will comment out everything about force. I will put a sign /**/ where I commented out
	//nvalues = 11;
	//nvalues = 24;	/**/
/////////////////////////////////DK: CHECK HERE. THIS IS PROBABLY THE REASON FOR THE BUG////////////////////////////////////////////////
	nvalues = 10;
	// 0-2  : dQ/dx
	// 3    : fatom
	// 4    : fcoup
	// 5    : fbias
	// 6    : fLwall
	// 7    : fUwall
	// 8-10 : dZ/dx (for z-restraining potential)

	size_peratom_cols = nvalues;	//number of columns in per atom array dQ/dx 
	peratom_freq = 1;		//frequency per-atom data is available at
	//nvalues = 24;
	// perform initial allocation of atom-based array
  	// register with Atom class
	grow_arrays(atom->nmax);
	atom->add_callback(0);
	
	// zero the array since dump may access it on timestep 0
  	// zero the array since a variable may access it before first run
	for (int i = 0; i < atom->nlocal; i++){
		for (int m = 0; m < nvalues; m++){
			array[i][m] = 0.0;
		}
	}


//DK: CHECK HERE
	//-----------------------
	// INCLUDE CHECK FOR UNITS IN LAMMPS!!!
	//-----------------------

	// Set up units
	// LAMMPS units wrt kj/mol - nm - ps
	// Set up units

  	/*if (force->boltz == 1.0){
		// LAMMPS units lj
    		p->cmd("setNaturalUnits");
  	} else {
    		double energyUnits=1.0;
    		double lengthUnits=1.0;
    		double timeUnits=1.0;
    		if (force->boltz == 0.0019872067){
			// LAMMPS units real :: kcal/mol; angstrom; fs
      			energyUnits=4.184;
      			lengthUnits=0.1;
      			timeUnits=0.001;
    		} else if (force->boltz == 8.617343e-5){
			// LAMMPS units metal :: eV; angstrom; ps
      			energyUnits=96.48530749925792;
      			lengthUnits=0.1;
      			timeUnits=1.0;
    		} else if (force->boltz == 1.3806504e-23){
			// LAMMPS units si :: Joule, m; s
      			energyUnits=0.001;
      			lengthUnits=1.e-9;
      			timeUnits=1.e-12;
    		} else if (force->boltz == 1.3806504e-16){
			// LAMMPS units cgs :: erg; cms;, s
      			energyUnits=6.0221418e13;
      			lengthUnits=1.e-7;
      			timeUnits=1.e-12;
    		} else if (force->boltz == 3.16681534e-6){
			// LAMMPS units electron :: Hartree, bohr, fs
      			energyUnits=2625.5257;
      			lengthUnits=0.052917725;
      			timeUnits=0.001;
    		} else error->all(FLERR,"Odd LAMMPS units, plumed cannot work with that");
    		p->cmd("setMDEnergyUnits",&energyUnits);
    		p->cmd("setMDLengthUnits",&lengthUnits);
    		p->cmd("setMDTimeUnits",&timeUnits);
  	}*/

	//currently only for metal, check what needs to be adjusted for this to work in other units!
//DK:MAKE SURE I CALCULATE HERE 
	/*if(force->boltz != 8.617343e-5){
		error->all(FLERR,"FOR FIX DAFED: proper units only checked for option metal, don't use other units at this time");
	}
	*/


	// Read fix parameters:
	//HERE:  EACH PROC READS ALL INPUT FILES, BUT ONLY PROC-0 WRITES STUFF
	// except for derivatives, there we have a different file for each proc!
  	int next=0;
	string dafedfile,cvfile,dcvfile;
	int cv_stride, dcv_stride;
	int ainmol;

  	for(int i=3;i<narg;++i){
    		if(!strcmp(arg[i],"outfile")) next=1;
    		else if(next==1){
      			/*if(universe->existflag == 1){
        			// Each replica writes an independent log file
        			//  with suffix equal to the replica id
        			char str_num[32], logFile[1024];
        			sprintf(str_num,".%d",universe->iworld);
        			strncpy(logFile,arg[i],1024-32);
        			strcat(logFile,str_num);
        			p->cmd("setLogFile",logFile);
        			next=0;
      			} else { */
        			// partition option not used
        			//p->cmd("setLogFile",arg[i]);
			if(me==0){
				cout << "The name of the output file is: "<<arg[i] << endl<<endl;
				string outfile = arg[i];
				of.open(outfile.c_str(),ofstream::out | ofstream::trunc);
				if (of.is_open()){
					of << "This is the DAFED output file!"<<endl;
					of << "This is checking narg: "<<narg << endl;
					of << "THis is checking what is arg[7]: " << arg[7] << endl;
				}
			}
        		next=0;
      			//}
    		}
    		else if(!strcmp(arg[i],"dafedfile"))next=2;
    		else if(next==2){
      			//p->cmd("setPlumedDat",arg[i]);
			cout <<me<< ": The name of the input file is: "<<arg[i] << endl<<endl;
			dafedfile = arg[i];
      			next=0;
    		}
		else if(!strcmp(arg[i],"totalsymf"))next=3;
    		else if(next==3){
  			parameter.nsfg = atoi(arg[i]);	
			cout << "number of symmetry function = " << parameter.nsfg << endl;
			next = 0;
		}
    		else error->all(FLERR,"syntax error in fix dafed - use 'fix name dafed dafedfile dafed.dat outfile dafed.out' ");
  	}
  	if(next==1) error->all(FLERR,"missing argument for outfile option");
  	if(next==2) error->all(FLERR,"missing argument for dafedfile option");
  	if(next==3) error->all(FLERR,"missing argument for totalsymf option");

	//get number of atoms
	int natoms=int(atom->natoms);	//this is the total number of atoms in the simulation

	
	//I added for the general use. 
	cout <<"parameter.nsfg" << parameter.nsfg << endl;
	parameter.RskappaLst = new double[parameter.nsfg];
	parameter.etaLst = new double[parameter.nsfg];
	cout <<"After para dynamic array" << endl;
	
	//read dafed input file
	DAFED::read_dafed_input(dafedfile,parameter,exvar,cvfile,cv_stride,dcvfile,dcv_stride);
	parameter.nmol = natoms/parameter.natm;	//no. of total mol in the system


	cout <<"parameter.natm " << parameter.natm << endl;
	cout <<"parameter.nsfg2CO " << parameter.nsfg2CO << endl;
	cout <<"parameter.nsfg2NN " << parameter.nsfg2NN << endl;
	cout <<"parameter.nnsfg3CO " << parameter.nsfg3CO << endl;
	cout <<"parameter.nsfg3NN " << parameter.nsfg3NN << endl;
	cout <<"parameter.nsfg2point " << parameter.nsfg2point << endl;
	cout <<"parameter.nsfg3point " << parameter.nsfg3point << endl;
	cout <<"parameter.nmol " << parameter.nmol << endl;
	cout <<"parameter.nnout " << parameter.nnout << endl;
	cout <<"parameter.center " << parameter.center << endl;
	cout <<"parameter.COvectype[0] " << parameter.COvectype[0] << endl;
	cout <<"parameter.COvectype[1] " << parameter.COvectype[1] << endl;
	cout <<"parameter.NNvectype[0] " << parameter.NNvectype[0] << endl;
	cout <<"parameter.NNvectype[1] " << parameter.NNvectype[1] << endl;
/*	
	for (int sy = 0; sy < parameter.nsfg;sy++)
	{
		cout << "Rskappa = " << parameter.RskappaLst[sy] << endl;
		cout << "etaLst = " << parameter.etaLst[sy] << endl;	
	}
*/

	//check if NN is used for structure determination
	if(parameter.useNN<0) error->all(FLERR,"!!Must define use_NN variable in dAFED input file!!");
	if(parameter.useNN==1){
		//read weights for the NN
		DAFED::read_weight(NNvar);
		//cout<<"proc "<<me<<" NNvar.input = "<<NNvar.input<<endl;
	}

	//seed for random number of extended variable (for initial velocities and langevin)
	int seed;
	seed = parameter.ex_seed;
	if(me==0){
		of<<endl<<"Seed for random number generator: "<<seed
			<<"  SAME ON EACH PROC, SHOULD BE OK SINCE INTEGRATION IS ONLY DONE BY PROC 0!!"<<endl;
	}
	DAFED::ran3(seed);
	//print some stuff to files
	if(me==0){
		of<<me<<"  First random number: "<<DAFED::ran3()<<endl;
		of<<endl;
		
		of<<"Name of the colvar file: "<<cvfile<<endl;
		init_colvar_file(cvfile);		//open and print header for colvar file

		if(parameter.useNN==1){
			init_nnout_file();		//open and print header for nnout file 
		}
		if(parameter.useMetadyn < 2){		//dafed and ufed
			init_dafedene_file();		//open and print header for dafed_ene file
		}

		//if(parameter.useMetadyn > 0){		//ufed and metadynamics
		init_bias_files();		//open files to record forces and bias potential on extended variable in UFED
		//}

		of<<"CV printed every "<<cv_stride<<" MD steps"<<endl;
	}

	nprint_cv = cv_stride;
	
	if(me==0){
		of<<endl;
		of<<"Name of the derivative file: "<<dcvfile<<endl;
	}
	//open derivatives files, one for each proc, and print headers
	init_derivatives_files(dcvfile);

	if(me==0){
		of<<"Derivatives printed every "<<dcv_stride<<" MD steps"<<endl;
	}
	nprint_deriv = dcv_stride;
	if(me==0){
        	of<<endl;	
		if(parameter.choose_exvar==0){
			if(parameter.nex==1){
				of<<"Using Q-A15 as collective variable"<<endl<<endl;
			}
			else if(parameter.nex==2){
				of<<"Using Q-BCC and Q-A15 as collective variables"<<endl<<endl;
			}
		}
		else if(parameter.choose_exvar==1){
			of<<"Using f(QBCC,QA15) as 1D collective variable"<<endl<<endl;
		}
		else if(parameter.choose_exvar==2){
			of<<"Using Q6 as collective variable"<<endl<<endl;
		}

		if(parameter.useMetadyn<2){
			of<<"Properties of the extended variables"<<endl;
			of<<"Number of extended variables: "<<parameter.nex<<endl;
		}
		else if(parameter.useMetadyn==2){
			if(parameter.lwall_cv != 1) error->all(FLERR,"With metadynamics bias can only be on collective variable, check para.lwall_cv!");
			for(int iex=0;iex<parameter.nex;iex++){
				if(exvar[iex].lwall_cv != 1) error->all(FLERR,"With metadynamics bias can only be on collective variable, check lwall_cv!");
				if(exvar[iex].uwall_cv != 1) error->all(FLERR,"With metadynamics bias can only be on collective variable, check uwall_cv!");
			}

			of<<"Properties of the collective variables -- Metadynamics simulation!!"<<endl;
			of<<"Number of collective variables: "<<parameter.nex<<endl;
		}
		else{
			error->all(FLERR,"useMetadyn has to be 0=dAFED, 1=ufed, or 2=metadynamics!!");
		}
	}
	int iex;
	double dt = update->dt;
	for(iex=0;iex<parameter.nex;iex++){
		//Convert to 
//		if (force->boltz == 0.0019872067)
//		{	// LAMMPS units real :: kcal/mol; angstrom; fs
//	    			exvar[iex].kappa = exvar[iex].kappa * 23.0609;	// eV-> kcal/mol
//  				exvar[iex].tau = exvar[iex].tau * 1000;		//fs
//				exvar[iex].m = exvar[iex].m * 23.0609*1000*1000;//ev*ps^2 -> kCal/mol*fs^2
//		}
		if(parameter.useMetadyn<2){
			if(me==0){
				if (force->boltz == 0.0019872067)
				{	// LAMMPS units real :: kcal/mol; angstrom; fs
					of<<"Kappa "<<iex<<": "<<exvar[iex].kappa<<" kCal/mol"<<endl;
					of<<"Tau "<<iex<<":   "<<exvar[iex].tau<<" fs"<<endl;
					//of<<"Gamma "<<iex<<": "<<exvar[iex].gamma<<" ps^-1"<<endl;
					of<<"T "<<iex<<":     "<<exvar[iex].temp<<" K"<<endl;
					of<<"Mass:   "<<exvar[iex].m<<" kCal/mol*fs^2"<<endl;
					of<<endl;
				}
				else
				{	// LAMMPS units real :: kcal/mol; angstrom; fs
					of<<"Kappa "<<iex<<": "<<exvar[iex].kappa<<" eV"<<endl;
					of<<"Tau "<<iex<<":   "<<exvar[iex].tau<<" ps"<<endl;
					//of<<"Gamma "<<iex<<": "<<exvar[iex].gamma<<" ps^-1"<<endl;
					of<<"T "<<iex<<":     "<<exvar[iex].temp<<" K"<<endl;
					of<<"Mass:   "<<exvar[iex].m<<" eV*ps^2"<<endl;
					of<<endl;
				}
			}
			//exvar[iex].m = exvar[iex].kappa*exvar[iex].tau*exvar[iex].tau/(4.0*M_PI*M_PI);
			//exvar[iex].m = exvar[iex].m*23.061;
			//initialize dafed
			exvar[iex].dt_respa = dt;  // In the current implementation, RESPA is done in the MD looop
			exvar[iex].thermo_work = 0.0;
			exvar[iex].dafed_work = 0.0;
			exvar[iex].dWold = 0.0;

			//set langevin parameter
			//exvar[iex].c1_lgvn = exp(-0.5*exvar[iex].gamma*dt);
			//exvar[iex].c2_lgvn = sqrt((1.0-exvar[iex].c1_lgvn*exvar[iex].c1_lgvn)*force->boltz*exvar[iex].temp/exvar[iex].m); 

			//initialize GGMT parameters
			//assumingt T and tau and n_respa are given
			//DK editted. to get eV from kCal/mol
			exvar[iex].ggmt.kT0 = force->boltz*exvar[iex].temp;
			double kT0 = exvar[iex].ggmt.kT0;
			double tau = exvar[iex].tau;
			exvar[iex].ggmt.Q1 = kT0*tau*tau;
			exvar[iex].ggmt.Q2 = 8.0/3.0*kT0*kT0*kT0*tau*tau;
			exvar[iex].ggmt.dt[0] = exvar[iex].dt_respa/((double)exvar[iex].ggmt.n_respa_ggmt);
			// Weights for the Suzuki-Yoshida decomposition:
			// These correspond to hte delta t_j from Liu et al. (2000)
			exvar[iex].ggmt.dt[1] = (1.0 / (2.0 - pow(2.0,(1.0/3.0)))) * exvar[iex].ggmt.dt[0];
			exvar[iex].ggmt.dt[2] = exvar[iex].ggmt.dt[0] - 2.0*exvar[iex].ggmt.dt[1];
			exvar[iex].ggmt.dt[3] = exvar[iex].ggmt.dt[1];

			exvar[iex].ggmt.eta1 = 1.0;
			exvar[iex].ggmt.eta2 = 1.0;
			exvar[iex].ggmt.v1 = -sqrt(kT0 / exvar[iex].ggmt.Q1);
			exvar[iex].ggmt.v2 = -sqrt(kT0 / exvar[iex].ggmt.Q2);

			if(me==0){
				of<<"GGMT "<<iex<<" - Q1:  "<<exvar[iex].ggmt.Q1<<endl;
				of<<"GGMT "<<iex<<" - Q2:  "<<exvar[iex].ggmt.Q2<<endl;
				of<<"GGMT "<<iex<<" - dt0:  "<<exvar[iex].ggmt.dt[0]<<endl;
				of<<"GGMT "<<iex<<" - dt1:  "<<exvar[iex].ggmt.dt[1]<<endl;
				of<<"GGMT "<<iex<<" - dt2:  "<<exvar[iex].ggmt.dt[2]<<endl;
				of<<"GGMT "<<iex<<" - dt3:  "<<exvar[iex].ggmt.dt[3]<<endl;
				of<<"GGMT "<<iex<<" - n_respa: "<<exvar[iex].ggmt.n_respa_ggmt;
				of<<endl;
			}
		}
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				exvar[iex].dvirial[i][j] = 0.0;
			}
		}
		exvar[iex].frestrain = 0.0;
		exvar[iex].vrestrain = 0.0;
	}


	if(me==0){
		of<<endl;
		if(parameter.lwall_cv == 0){
			of<<"Lower wall on SUM OF EXTENDED VARIABLES"<<endl;
		}
		else if(parameter.lwall_cv == 1){
			of<<"Lower wall on SUM OF COLLECTIVE VARIABLES"<<endl;
		}
		else error->all(FLERR,"set value of lwall_cv to apply lower wall to sum of either 0=exvar or 1=cv");
		of<<"Lower wall value for the sum of extended variables, slimit = "<<parameter.lwall<<endl;
		of<<"Lower wall potential parameters V(s)=k*((s-slimit)/eps)^n:"<<endl;
		of<<"n = "<<parameter.lwall_n<<"  eps = "<<parameter.lwall_eps<<"  k = "<<parameter.lwall_k<<endl;
		of<<endl;
		for(iex=0;iex<parameter.nex;iex++){
			if(exvar[iex].lwall_cv == 0){
				of<<"Lower wall constraint applied to extended variable (exvar)"<<endl;
			}
			else if(exvar[iex].lwall_cv == 1){
				of<<"Lower wall constraint applied to collective variable (cv)"<<endl;
			}
			else error->all(FLERR,"set value of ex_lwall_cv to apply lower wall constraint to either 0=exvar or 1=cv");
			of<<"Lower wall value for extended variables "<<iex<<", slimit = "<<exvar[iex].lwall<<endl;
			of<<"Lower wall potential parameters V(s)=k*((s-slimit)/eps)^n:"<<endl;
			of<<"n = "<<exvar[iex].lwall_n<<"  eps = "<<exvar[iex].lwall_eps<<"  k = "<<exvar[iex].lwall_k<<endl;
			of<<endl;
			exvar[iex].lwall_pref = -exvar[iex].lwall_k*exvar[iex].lwall_n;
			for(int i=0;i<exvar[iex].lwall_n;i++){
				exvar[iex].lwall_pref /= exvar[iex].lwall_eps;
			}


			if(exvar[iex].uwall_cv == 0){
				of<<"Upper wall constraint applied to extended variable (exvar)"<<endl;
			}
			else if(exvar[iex].uwall_cv == 1){
				of<<"Upper wall constraint applied to collective variable (cv)"<<endl;
			}
			else error->all(FLERR,"set value of ex_uwall_cv to apply lower wall constraint to either 0=exvar or 1=cv");
			of<<"Upper wall value for extended variables "<<iex<<", slimit = "<<exvar[iex].uwall<<endl;
			of<<"Upper wall potential parameters V(s)=k*((s-slimit)/eps)^n:"<<endl;
			of<<"n = "<<exvar[iex].uwall_n<<"  eps = "<<exvar[iex].uwall_eps<<"  k = "<<exvar[iex].uwall_k<<endl;
			of<<endl;
			exvar[iex].uwall_pref = -exvar[iex].uwall_k*exvar[iex].uwall_n;
			for(int i=0;i<exvar[iex].uwall_n;i++){
				exvar[iex].uwall_pref /= exvar[iex].uwall_eps;
			}
		}
		if(parameter.use_restraint==1){
			if(parameter.choose_exvar != 1){
				//error->all(FLERR,"Additional restraint can only be used together with choose_exvar=1!!");
				cerr<<"\n\n\t!!ERROR!! Additional restraint can only be used together with choose_exvar=1!!"<<endl;
				cerr<<"\tExiting programme..."<<endl<<endl;;
				exit(1);
			}
			of<<endl;
			of<<"Using additional restraining potential on collective variable S=BCC+A15 (tube)"<<endl;
			of<<"Restraining value, slimit = "<<parameter.restraint_value<<endl;
			of<<"Restraining potential parameters V(s)=k*((s-slimit)/eps)^n:"<<endl;
			of<<"n = "<<parameter.restraint_n<<"  eps = "<<parameter.restraint_eps<<"  k = "<<parameter.restraint_k<<endl;
			of<<endl;
			
		}
		of<<endl;
	}
	if(parameter.use_restraint==1){
		parameter.restraint_pref = -parameter.restraint_k * parameter.restraint_n;
		for(int i=0; i<parameter.restraint_n; i++){
			parameter.restraint_pref /= parameter.restraint_eps;
		}
	}


	//set langevin parameter
	//double gamma = exprop[2];
	//double dt = update->dt;
	//c1_lgvn = exp(-0.5*gamma*dt);
	//c2_lgvn = sqrt((1.0-c1_lgvn*c1_lgvn)*force->boltz*T_s/m_s); 


	//this sets that the fix contributes to the virial and that it has a scalar value
	//virial_flag=1;		// 1 if the fix contributes to the virial!  HERE WE DO NOT CORRECTLY DO THIS YET!!
	virial_global_flag = virial_peratom_flag = 1;
	//scalar_flag = 1;

	// Define compute to calculate potential energy
	char *id_pe = (char *) "thermo_pe";
	int ipe = modify->find_compute(id_pe);
	c_pe = modify->compute[ipe];
	// Trigger computation of potential energy every step
	c_pe->addstep(update->ntimestep+1);

	//set parameters for the cutoff and symmetry functions
	//this is the larger cutoff for the symmetry functions
//	parameter.rmin0 = 9.8;	
//	parameter.rmax0 = 10;
	//Need to fix this part.
  	//parameter.rmin0 = cutoff_user-0.2;	
  	//parameter.rmax0 = cutoff_user;
	//this is the smaller cutoff for the Steinhardt parameters
	//parameter.rmin1 = 3.8;
	//parameter.rmax1 = 4.0;
	
	parameter.npairs = 0;
	








	if(parameter.useNN==1){
		//for the symmetry function - needs to be read in from file
	/*	parameter.nsfg=24;
		parameter.nsfg2CO=4;	//no. of G3 symmetry functions
		parameter.nsfg2NN=4;	//no. of G3 symmetry functions
		parameter.nsfg3CO=4;	//no. of G3 symmetry functions
		parameter.nsfg3NN=4;	//no. of G3 symmetry functions
		parameter.nsfg2point=4;	//no. of G3 symmetry functions
		parameter.nsfg3point=4;	//no. of G3 symmetry functions
		parameter.nmol = natoms/parameter.natm;	//no. of total mol in the system
		parameter.center = 3;
		parameter.nnout = 6;	
		parameter.COvectype[0] = 3;
		parameter.COvectype[1] = 2;
		parameter.NNvectype[0] = 4;
		parameter.NNvectype[1] = 5;
	*/	if(parameter.nsfg != NNvar.input){
			error->all(FLERR,"Number of calculated Steinhard parameter and symmetry functions NOT THE SAME as NN input values!");
		}

		if(me==0){
			int isym;
			of<<"Total number of input functions for NN: "<<parameter.nsfg<<endl;
			//of<<"Using "<<parameter.nstein<<" Steinhardt parameters: q6, q7, q8"<<endl;
			of<<"Values of the symmetry function parameters"<<endl;
			of<<"sfg   Rskappa   eta"<<endl;
			for(isym=0;isym<parameter.nsfg;isym++){
				of<<isym<<"    "<<Rskappa[isym]<<"    "<<eta[isym]<<endl;
			}
			of<<endl;
		}
	}
	cout << "ERROR CHECKER" << endl;
	// max # of owned+ghost in arrays on this proc	
	nmax = atom->nmax;
	cout << "nmax = " << nmax << endl;
	//allocate memory for extended variable derivatives
	//DK: Make sure I did not overwrite this ev with myexvar_temp 
	for(iex=0;iex<parameter.nex;iex++){
		//allocate derivatives
		exvar[iex].dQ = new double* [nmax];
		for(int j=0;j<nmax;j++){
			exvar[iex].dQ[j] = new double [3];
		}
	}

	//allocate memory for atoms data class
	//initialize molecules data structure
	molecules = new DAFED::CAtom[nmax];
	molec = new DAFED::Cmolpoint[nmax];
	if(parameter.useNN == 1){
		for(int i=0;i<nmax;i++){
			molecules[i].sfg = new double[parameter.nsfg];
			molecules[i].NNout = new double[NNvar.output];
			molecules[i].NNgrad = new double[NNvar.output*NNvar.input];
			molec[i].atom_ID = new int[parameter.natm];
		}
	}
	//stuff from fcc/orient
  	nbr = (Nbr *) memory->smalloc(nmax*sizeof(Nbr),"fix/dafed:nbr");

	//if(parameter.useNN == 1){
	//	comm_forward = 166 + (NNvar.output*NNvar.input);	//no of values that is being communicated for each atoms!
	//}
	//else{
	//	comm_forward = 166;
	//}
	memory->create(NNout,nmax,parameter.nnout,"dafedgen:NNout");
	
	
	
	
	array_atom = NNout;

	comm_forward = 6 + (NNvar.output*NNvar.input);
	//comm_reverse = 3;
	//comm_reverse = 9;
	comm_reverse = 3+parameter.nnout;

	//initialize values of virial, in first step there is no contribution from exvar, since exvar.f=0
	for(int i=0;i<6;i++){
		Fix::virial[i] = 0.0;
	}

	//check if additional bias potential is used and setup corresponding histogram
	if(parameter.useMetadyn > 0){
		if(me == 0){
			setup_ufedbias();
			count_updatebias = 0;
		}
		count_printbias = 0;
		nprint_bias = parameter.bias_stride;

		if(me == 0){
			if(parameter.bias_read == 1){		//read in existing bias potential
				of<<"Reading in bias potential from file: "<<parameter.bias_readfile;
				of<<endl<<endl;
				read_ufedbias();
			}
		}

	}

	// begin debug
	//ostringstream ssdebug;
	//string debugfile;
	//ssdebug<<me;
	//debugfile = "debug_virial-"+ssdebug.str();
	//of_debug_virial.open(debugfile.c_str(),ofstream::out | ofstream::trunc);
	// end debug

	//something for constant force
	biasforce_sign = 1.0;

	MPI_Barrier(world);



}


// ---------------------------------------------------------------------- 
// destructor
// ---------------------------------------------------------------------- 
FixDafedgen::~FixDafedgen()
{
	int iex,i,j;
	//free memory of derivatives of extended variables
	memory->destroy(NNout);

	delete [] parameter.RskappaLst;
	delete [] parameter.etaLst;
	for(iex=0;iex<parameter.nex;iex++){
		for(j=0;j<nmax;j++){
			delete [] exvar[iex].dQ[j];
		}
		delete [] exvar[iex].dQ;
	}

	if(parameter.useNN == 1){
		for(i=0;i<nmax;i++){
			delete [] molecules[i].sfg;
			delete [] molecules[i].NNout;
			delete [] molecules[i].NNgrad;
			delete [] molec[i].atom_ID;
		}
	}
	delete [] molecules;
	delete [] molec;
	if(parameter.useMetadyn > 0){
		if(me == 0){
			for(i=0; i<Nlattot; i++){
				delete [] bhist[i].x;
			}
			delete [] bhist;
			for(i=0;i<Nvertex;i++){
				delete [] nblat[i].diff;
			}
			delete [] nblat;
		}
	}

	memory->sfree(nbr);

	atom->delete_callback(id,0);
	memory->destroy(array);

}


// ---------------------------------------------------------------------- 
//determines when the fix is called during the timestep (required)
// ---------------------------------------------------------------------- 
int FixDafedgen::setmask()
{
  // set with a bitmask how and when apply the force from dafed
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

// ---------------------------------------------------------------------- 
//initialization before a run (optional)
// ---------------------------------------------------------------------- 
void FixDafedgen::init()
{
	//cout<<"In the init of fix_dafed"<<endl;
	if (strcmp(update->integrate_style,"respa") == 0)
    		nlevels_respa = ((Respa *) update->integrate)->nlevels;

	// need a full neighbor list
  	// perpetual list, built whenever re-neighboring occurs
	//OLD code here. no longer runable (6 lines from here. 7th line is a replacement)
  	//int irequest = neighbor->request(this,instance_me);
  	//neighbor->requests[irequest]->pair = 0;
  	//neighbor->requests[irequest]->fix = 1;
  	//neighbor->requests[irequest]->half = 0;
  	//neighbor->requests[irequest]->full = 1;
	//neighbor->requests[irequest]->ghost = 1;	//if you need a neighbour list for the ghost atoms
	neighbor->add_request(this, NeighConst::REQ_FULL);


}

// ---------------------------------------------------------------------- 
// init list needed for communciation
// ---------------------------------------------------------------------- 
void FixDafedgen::init_list(int id, NeighList *ptr)
{
  list = ptr;
}


// ---------------------------------------------------------------------- 
//called immediately before the 1st timestep and after forces are computed (optional)
// ------------------------------------------------------------------------
void FixDafedgen::setup(int vflag)
{
	
	//cout<<"In setup of fix_dafed"<<endl;
	if (strcmp(update->integrate_style,"verlet") == 0){
		post_force(vflag);   //so, here it calls the post_force function, BEFORE the first time step!!
		end_of_step();       //and to the end_of_step function, BEFORE the first time step!!
	}

	else {
		((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
		post_force_respa(vflag,nlevels_respa-1,0);
		((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
		end_of_step();       //and to the end_of_step function, BEFORE the first time step!!

	}

}



// ------------------------------------------------------------------------
// called after pair & molecular forces are computed and communicated (optional)
// THIS IS THE MAIN FUNCTION HERE! 
// ------------------------------------------------------------------------
void FixDafedgen::post_force(int vflag)
{
	int natoms=int(atom->natoms);	//total number of atoms
	int nlocal=int(atom->nlocal);	//total number of atoms on this proc
	int i,j;
	int step=update->ntimestep;
	double dt = update->dt;
	double time;
	int nex=parameter.nex;
	int iex;
	//double vex_old[2], xex_old[2];	//store old positions and velocities if you run into the lower boundary
	double kt;
	double exvarsum;		//determine sum of either exvar or cv
	double pos_constraint;		//position used to determined force from constraint, either from exvar or cv

	double **x = atom->x;
	double **f = atom->f;
	tagint *tag = atom->tag;

	//these are the parameters for the wall potential Vwall = k((s-slimit)/eps)^n
	//this is a lower wall for \sum s_i (sum over extended variables)
	double fwall=0.0;
	double epswall=parameter.lwall_eps;
	double kwall=parameter.lwall_k;
	int nwall=parameter.lwall_n;
	double prefwall;

	prefwall = -kwall*nwall;
	for(i=0;i<nwall;i++){
		prefwall /= epswall;
	}


	//this is being used in the function that calculates value of collective variable
	lmpneigh.ntot = atom->natoms;
	lmpneigh.nall = atom->nlocal + atom->nghost;
	
	if (lmpneigh.nall > nmax) {	// re-allocate memory if no. of atoms is larger than nmax!

		//first delete the variables
		memory->destroy(NNout);
		for(iex=0;iex<parameter.nex;iex++){
			for(j=0;j<nmax;j++){
				delete [] exvar[iex].dQ[j];
			}
			delete [] exvar[iex].dQ;
		}


		if(parameter.useNN == 1){
			for(i=0;i<nmax;i++){
				delete [] molecules[i].sfg;
				delete [] molecules[i].NNout;
				delete [] molecules[i].NNgrad;
				delete [] molec[i].atom_ID;
			}
		}
		delete [] molecules;
		delete [] molec;
		//new size considering all the atoms on this proc
    		nmax = lmpneigh.nall;
		for(iex=0;iex<parameter.nex;iex++){
			//allocate derivatives
			exvar[iex].dQ = new double* [nmax];
			for(int j=0;j<nmax;j++){
				exvar[iex].dQ[j] = new double [3];
			}
		}

		//initialize molecules data structure
		molecules = new DAFED::CAtom[nmax];
		molec = new DAFED::Cmolpoint[nmax];
		memory->create(NNout,nmax,parameter.nnout,"dafedgen:NNout");
		array_atom = NNout;
		if(parameter.useNN == 1){
			for(i=0;i<nmax;i++){
				molecules[i].sfg = new double[parameter.nsfg];
				molecules[i].NNout = new double[NNvar.output];
				molecules[i].NNgrad = new double[NNvar.output*NNvar.input];
				molec[i].atom_ID = new int[parameter.natm];
			}
		}

		memory->destroy(nbr);
    		nbr = (Nbr *) memory->smalloc(nmax*sizeof(Nbr),"orient/fcc:nbr");

	}

	for(iex=0;iex<nex;iex++){
		exvar[iex].Q = 0.0;		//this is a 'global' value that is the same on each proc
	}

	
	//get Q6 value and derivatives (this is AFTER atoms have been moved, so Q(r(t+dt))
	//except when this is called from setup(), for initial configuration
	
	//DAFED::get_aq6(atom->x,natoms,box,dQ6_pairs,Q6_pairs);
	//direct output of the NN
	//DAFED::get_Q_short(atom->x,natoms,box,exvar,parameter,sf2Rs,sf2eta,sf3kappa,NNvar,Qall);
	//normalised output of the NN
	//DAFED::get_Qnorm_short(atom->x,natoms,box,exvar,parameter,sf2Rs,sf2eta,sf3kappa,NNvar,Qall);
	//using Steinhardt parameters and symmetry functions as input for NN
	//DAFED::get_Q_short_3(atom->x,natoms,box,exvar,parameter,sf2Rs,sf2eta,sf3kappa,NNvar,Qall);
	//using only the 3 Steinhardt parameters
	//DAFED::get_Q_short_stein(atom->x,natoms,box,exvar,parameter,NNvar,Qall);
	
	//for the parallel version
	// get the values of the extended variables and their derivatives!
	if(parameter.choose_exvar == 0){
		cout << "ERROR"; //get_exvar_values();
	}
	else if(parameter.choose_exvar == 1){
		//get_pathcv_values_test();
		get_pathcv_values_gen();
		//get_pathcv_values_Xpoints();
	}	       
	else if(parameter.choose_exvar == 2){
		//get_Q6_values();
		cout << "ERROR"; //get_exvar_values();
	}
	else error->all(FLERR,"In post_force: this is weird, variable choose_exvar can only be 0, 1 or 2!");
	MPI_Barrier(world);	//not sure if this is needed since there is an MPI_Allreduce in get_exvar_values
	//exit(1);	
//	get_numerical_derivatives();
//	exit(1);
	//before the first integration step
	if(step==0){

		//initialize virial
		for(i=0;i<6;i++){
			Fix::virial[i] = 0.0;
		}

		exvarsum = 0.0;
		for(iex=0;iex<nex;iex++){
			if(me==0){			//only modify extended variable on root proc
				//set initial value of extended variable
				exvar[iex].x = exvar[iex].Q;
				if(parameter.lwall_cv == 0){
					exvarsum += exvar[iex].x;
				}
				else if(parameter.lwall_cv == 1){
					exvarsum += exvar[iex].Q;
				}
				else error->all(FLERR,"para.lwall_cv not properly set!");
				if(parameter.useMetadyn < 2){	//dafed and ufed
					//set initial velocity
					exvar[iex].v = DAFED::gauss()*sqrt(force->boltz*exvar[iex].temp/exvar[iex].m);  // [1/fs]
					//inital force
					exvar[iex].fcoup = exvar[iex].kappa*(exvar[iex].Q - exvar[iex].x); //if we set these equal, the inital force is 0
					exvar[iex].f = exvar[iex].fcoup;
					if(parameter.useMetadyn > 0){
						exvar[iex].f += exvar[iex].fbias;
					}
					exvar[iex].fatom = -exvar[iex].fcoup;
				}
				//inital force of lower wall is 0.0 since right now we only allow to start simulatio within boundaries
				exvar[iex].fwall = 0.0;
				exvar[iex].fLwall = 0.0;
				exvar[iex].fUwall = 0.0;
				if(parameter.useMetadyn == 2){	//metdaynamics
					exvar[iex].fatom = exvar[iex].fbias;
					exvar[iex].fatom += exvar[iex].fwall;
					exvar[iex].fatom += exvar[iex].fLwall;
					exvar[iex].fatom += exvar[iex].fUwall;

					exvar[iex].f = 0.0;
					exvar[iex].v = 0.0;
					exvar[iex].fcoup = 0.0;
				}

			}
			MPI_Barrier(world);	
			
			//communicate position, veolcity and force on extended variable
			MPI_Bcast(&exvar[iex].x,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].v,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].f,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fcoup,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fbias,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fwall,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fLwall,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fUwall,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fatom,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvarsum,1,MPI_DOUBLE,0,world);
			
			//add to forces of all the atoms
			//in UFED this is only the part from the coupling
			for(i=0;i<nlocal;i++){
				for(j=0;j<3;j++){
					//atom->f[i][j] -= exvar[iex].fcoup*exvar[iex].dQ[i][j];
					atom->f[i][j] += exvar[iex].fatom*exvar[iex].dQ[i][j];		//this should be ok, -dV/dq*dq/dr
				}
			}
			
			//update virial
			//in UFED this is only the part from the coupling
			//Fix::virial[0] += -exvar[iex].fcoup*exvar[iex].dvirial[0][0];				
			//Fix::virial[1] += -exvar[iex].fcoup*exvar[iex].dvirial[1][1];				
			//Fix::virial[2] += -exvar[iex].fcoup*exvar[iex].dvirial[2][2];				
			//Fix::virial[3] += -exvar[iex].fcoup*exvar[iex].dvirial[0][1];				
			//Fix::virial[4] += -exvar[iex].fcoup*exvar[iex].dvirial[0][2];				
			//Fix::virial[5] += -exvar[iex].fcoup*exvar[iex].dvirial[1][2];

			Fix::virial[0] += exvar[iex].fatom*exvar[iex].dvirial[0][0];				
			Fix::virial[1] += exvar[iex].fatom*exvar[iex].dvirial[1][1];				
			Fix::virial[2] += exvar[iex].fatom*exvar[iex].dvirial[2][2];				
			Fix::virial[3] += exvar[iex].fatom*exvar[iex].dvirial[0][1];				
			Fix::virial[4] += exvar[iex].fatom*exvar[iex].dvirial[0][2];				
			Fix::virial[5] += exvar[iex].fatom*exvar[iex].dvirial[1][2];

		}	
		//begin debug
		//of_debug_virial<<setw(7)<<step<<"  ";
		//of_debug_virial<<setprecision(4)<<fixed<<time<<"  ";
		//for(i=0;i<6;i++){
		//	of_debug_virial<<setprecision(6)<<scientific<<Fix::virial[i]<<"  ";
		//}
		//of_debug_virial<<endl;
		//end debug
		
		if(exvarsum < parameter.lwall){	//this is for putting a hard wall at s1+s2 > lwall
			of<<"\n\n\t!!ERROR!! Sum of extended/collective variables below lower wall!\n";
			of<<"Lower wall = "<<parameter.lwall<<endl;
			of<<"Extended/collective variables:"<<endl;
			for(iex=0;iex<nex;iex++){
				of<<"s"<<iex<<" = "<<exvar[iex].x
					<<"cv"<<iex<<" = "<<exvar[iex].Q
					<<endl;
			}
			of<<"exvarsum = "<<exvarsum<<endl<<endl;
			error->all(FLERR,"values of extended/collective variables too small (below lower_wall), check log file and initial configuration! ");
		}
		for(iex=0;iex<nex;iex++){
			if(exvar[iex].lwall_cv == 0){
				pos_constraint = exvar[iex].x;
			}
			else if(exvar[iex].lwall_cv == 1){
				pos_constraint = exvar[iex].Q;
			}
			else error->all(FLERR,"lwall_cv not properly set!");
			//if(exvar[iex].x < exvar[iex].lwall){	//this is a lower wall for each exvar
			if(pos_constraint < exvar[iex].lwall){	//this is a lower wall for each exvar
				//of<<"\n\n\t!!ERROR!! Extended variables "<<iex<<" below lower wall!\n";
				of<<"\n\n\t!!WARNING!! Extended variables "<<iex<<" below lower wall!\n";
				of<<"Lower wall = "<<exvar[iex].lwall<<endl;
				of<<"Extended variable:"<<endl;
				of<<"s"<<iex<<" = "<<exvar[iex].x<<endl;
				of<<"Collective variable:"<<endl;
				of<<"cv"<<iex<<" = "<<exvar[iex].Q<<endl;
				//error->all(FLERR,"value of extended variable too small (below lower_wall), check log file and initial configuration! ");
			}

			if(exvar[iex].uwall_cv == 0){
				pos_constraint = exvar[iex].x;
			}
			else if(exvar[iex].uwall_cv == 1){
				pos_constraint = exvar[iex].Q;
			}
			else error->all(FLERR,"uwall_cv not properly set!");
			//if(exvar[iex].x > exvar[iex].uwall){	//this is a upper wall for each exvar
			if(pos_constraint > exvar[iex].uwall){	//this is a upper wall for each exvar
				of<<"\n\n\t!!WARNING!! Extended variables "<<iex<<" above upper wall!\n";
				of<<"Upper wall = "<<exvar[iex].uwall<<endl;
				of<<"Extended variable:"<<endl;
				of<<"s"<<iex<<" = "<<exvar[iex].x<<endl;
				of<<"Collective variable:"<<endl;
				of<<"cv"<<iex<<" = "<<exvar[iex].Q<<endl;
				//error->all(FLERR,"value of extended variable too large (above upper_wall), check log file and initial configuration! ");
			}
		}

		time = update->atime + (update->ntimestep-update->atimestep)*update->dt;
		//write some stuff if you are proc 0
		if(me==0){
			print_colvar(step,time);	//print exvar/colvar values to file
			
			if(parameter.useNN == 1){
				print_nnout(step,time);
			}

			if(parameter.useMetadyn < 2){		//dafed and ufed
				of_ene<<setw(7)<<step<<"  ";
				of_ene<<setprecision(4)<<fixed<<time<<"  ";
				for(iex=0;iex<nex;iex++){
					print_dafed(exvar[iex]);
				}
				of_ene<<endl;
			}
			//print all forces
			of_biasforce<<setw(7)<<step<<"  ";
			of_biasforce<<setprecision(4)<<fixed<<time<<"  ";
			for(iex=0;iex<nex;iex++){
				print_biasforce(exvar[iex]);
			}
			of_biasforce<<endl;

			if(parameter.useMetadyn > 0){		//ufed and metadynamics
				//print ufed bias potential
				print_biaspot();
			}
					
		
		}
		//wait until proc 0 has written all this
		MPI_Barrier(world);
		//Every proc writes derivatives wrt its atoms
		of_deriv<<"#Time step: "<<step<<endl;
		for(i=0;i<nlocal;i++){
			of_deriv<<setw(6)<<atom->tag[i];
			for(iex=0;iex<nex;iex++){
				for(j=0;j<3;j++){
					of_deriv<<setw(18)<<setprecision(6)<<scientific<<exvar[iex].dQ[i][j];
				}
			}
			of_deriv<<endl;
		}

		//get numerical derivatives and compare with analytical ones
		//get_numerical_derivatives();
		//exit(1);
	}
	else{
		if(me==0){		//update some stuff on the extended variable
			if(parameter.useMetadyn > 0){
				//update UFED bias potential
				if(count_updatebias == parameter.gauss_freq){
					update_ufedbias();
					count_updatebias = 0;
				}
				count_updatebias += 1;
			}

			exvarsum = 0.0;
			for(iex=0;iex<nex;iex++){
				if(parameter.useMetadyn < 2){	//integrate extended variable, dafed and ufed
					//Langevin integrator
					//first half langevin step
					//exvar[iex].v = exvar[iex].c1_lgvn*exvar[iex].v + exvar[iex].c2_lgvn*DAFED::gauss();
					//first half step of velocities with forces at time t
					//exvar[iex].v += 0.5*exvar[iex].f*dt/exvar[iex].m;
					//update positions
					//exvar[iex].x += exvar[iex].v*dt;

					//GGMT integrator
					//integrate from 0 to dt/2
					integrate_ggmt(exvar[iex]);
					//update velocities
					exvar[iex].v += 0.5*exvar[iex].dt_respa*exvar[iex].f/exvar[iex].m;
					//update positions
					exvar[iex].x += exvar[iex].dt_respa * exvar[iex].v;
					//sum positions for putting a wall on the sum
					if(parameter.lwall_cv == 0){		//wall on \sum exvar_i
						exvarsum += exvar[iex].x;
					}
					else if(parameter.lwall_cv == 1){	//wall on \sum cv_i
						exvarsum += exvar[iex].Q;
					}
					else error->all(FLERR,"para.lwall_cv not properly set!");
			
					// Integrate the thermostat work (for now only simple rectangular scheme)
					// in principle here we have kt and v1, v2 and dt/2
					kt = exvar[iex].v*exvar[iex].v*exvar[iex].m;
					exvar[iex].thermo_work -= 0.5*exvar[iex].dt_respa *
						kt*(exvar[iex].ggmt.v1 + exvar[iex].ggmt.v2*(exvar[iex].ggmt.kT0 + kt/3.0));

					//update the force (force at t+dt)
					exvar[iex].fcoup = exvar[iex].kappa*(exvar[iex].Q - exvar[iex].x);
					exvar[iex].f = exvar[iex].fcoup;
					exvar[iex].fatom = -exvar[iex].fcoup;
				}
				else if(parameter.useMetadyn == 2){	//metadynamics, x is value of cv
					exvar[iex].x = exvar[iex].Q;
					exvarsum += exvar[iex].x;
					exvar[iex].fatom = 0.0;
				}
				else{
					cerr<<"\n\n\t!!ERROR!! Something weird with useMetadyn parameter in post_force!!"<<endl;
					cerr<<"\tShould be 0=dafed, 1=ufed, or 2=metadynamics\n";
					cerr<<"\tparameter.useMetadyn = "<<parameter.useMetadyn<<endl;
					cerr<<"\tExiting programme...\n";
					exit(1);
				}


				//check for additional forces due to lower/upper wall
				//lower wall
				//check if bias is on exvar or cv
				exvar[iex].fLwall = 0.0;
				if(exvar[iex].lwall_cv == 0){
					pos_constraint = exvar[iex].x;
				}
				else if(exvar[iex].lwall_cv == 1){
					pos_constraint = exvar[iex].Q;
				}
				else error->all(FLERR,"lwall_cv not properly set!");	
				//if(exvar[iex].x < exvar[iex].lwall){
				if(pos_constraint < exvar[iex].lwall){
					double diffwall, diffwall_prod;
					//diffwall=exvar[iex].x-exvar[iex].lwall;
					diffwall=pos_constraint-exvar[iex].lwall;
					diffwall_prod = 1.0;
					for(int iwall=0; iwall<(exvar[iex].lwall_n-1);iwall++){
						diffwall_prod *= diffwall;
					}
					exvar[iex].fLwall = exvar[iex].lwall_pref*diffwall_prod;
					//if(parameter.useMetadyn < 2){		//dafed and ufed
					if(exvar[iex].lwall_cv == 0){		//constraint on exvar
						// add to total force of exvar
						exvar[iex].f += exvar[iex].fLwall;
					}
					//else if(parameter.useMetadyn == 2){	//metadynamics
					else if(exvar[iex].lwall_cv == 1){	//constraint on cv
						exvar[iex].fatom += exvar[iex].fLwall;
					}
					else{
						cerr<<"\n\n\t!!ERROR!! Lwall: lwall_cv = "<<exvar[iex].lwall_cv<<endl;
						cerr<<"\tExiting programme...\n";
						exit(1);
					}
					//JR: bias-wall start
					//also update the bias potential for this current position! Avoid trapping
					//if(parameter.useMetadyn > 0){
					//	//update bias potential
					//	update_ufedbias();
					//}
					//JR: bias-wall end

				}
				//upper wall
				//check if bias is on exvar or cv
				exvar[iex].fUwall = 0.0;
				if(exvar[iex].uwall_cv == 0){
					pos_constraint = exvar[iex].x;
				}
				else if(exvar[iex].uwall_cv == 1){
					pos_constraint = exvar[iex].Q;
				}
				else error->all(FLERR,"uwall_cv not properly set!");
				//if(exvar[iex].x > exvar[iex].uwall){
				if(pos_constraint > exvar[iex].uwall){
					double diffwall, diffwall_prod;
					//diffwall=exvar[iex].x-exvar[iex].uwall;
					diffwall=pos_constraint-exvar[iex].uwall;
					diffwall_prod = 1.0;
					for(int iwall=0; iwall<(exvar[iex].uwall_n-1);iwall++){
						diffwall_prod *= diffwall;
					}
					exvar[iex].fUwall = exvar[iex].uwall_pref*diffwall_prod;
					//if(parameter.useMetadyn < 2){		//dafed and ufed
					if(exvar[iex].uwall_cv == 0){		//constraint on exvar
						// add to total force
						exvar[iex].f += exvar[iex].fUwall;
					}
					//else if(parameter.useMetadyn == 2){	//metadynamics
					else if(exvar[iex].uwall_cv == 1){	//constraint on cv
						exvar[iex].fatom += exvar[iex].fUwall;
					}
					else{
						cerr<<"\n\n\t!!ERROR!! Uwall: uwall_cv = "<<exvar[iex].uwall_cv<<endl;
						cerr<<"\tExiting programme...\n";
						exit(1);
					}
					//JR: bias-wall start
					//also update the bias potential for this current position! Avoid trapping
					//if(parameter.useMetadyn > 0){
					//	//update bias potential
					//	update_ufedbias();
					//}
					//JR: bias-wall end

				}

			}
			//additional force from bias potential
			//if(parameter.useMetadyn > 0){		//ufed and metadyanmics
			//	get_bias_force();
			//}

			//now check if the sum is below lwall, exvarsum is either \sum exvar_i or \sum cv_i
			if(exvarsum < parameter.lwall){	
				
				//determine additional contribution to the force from the lower bound
				double diffwall, diffwall_prod;
				diffwall=exvarsum-parameter.lwall;
				diffwall_prod = 1.0;
				for(int iwall=0; iwall<(nwall-1);iwall++){
					diffwall_prod *= diffwall;
				}
				fwall = prefwall*diffwall_prod;
				for(iex=0;iex<nex;iex++){
					exvar[iex].fwall = fwall;
				}
				//JR: bias-wall start
				//also update the bias potential for this current position!  Avoid trapping
				//if(parameter.useMetadyn > 0){
				//	//update bias potential
				//	update_ufedbias();
				//}
				//JR: bias-wall end
			}
			else{
				for(iex=0;iex<nex;iex++){
					exvar[iex].fwall = 0.0;
				}
			}
			//additional force from bias potential
			if(parameter.useMetadyn > 0){		//ufed and metadyanmics
				get_bias_force();
				//get_bias_fixedforce(); 		//JR this is a quick fix for a constant force!
			}

			for(iex=0;iex<nex;iex++){
				if(parameter.useMetadyn < 2){		//dafed and ufed
					//follow 'normal' GGMT integration
					//additional forces from ufed bias potential
					if(parameter.useMetadyn > 0){
						exvar[iex].f += exvar[iex].fbias;
					}
					//additional force form lower wall
					if(parameter.lwall_cv == 0){			//on exvar
						exvar[iex].f += exvar[iex].fwall;
					}
					else if(parameter.lwall_cv == 1){		//on cv
						exvar[iex].fatom += exvar[iex].fwall;
					}
					else error->all(FLERR,"para.lwall_cv not properly set!");
			
					//integrate from dt/2 to dt
					exvar[iex].v += 0.5*exvar[iex].dt_respa*exvar[iex].f/exvar[iex].m;
					integrate_ggmt(exvar[iex]);

					// Integrate the thermostat work (for now only simple rectangular scheme)
					// in principle here we have kt and v1, v2 and dt
					kt = exvar[iex].v*exvar[iex].v*exvar[iex].m;
					exvar[iex].thermo_work -= 0.5*exvar[iex].dt_respa *
						kt*(exvar[iex].ggmt.v1 + exvar[iex].ggmt.v2*(exvar[iex].ggmt.kT0 + kt/3.0));

					// Integrate DAFED work (triangle scheme)
					// only work from coupling force
					double dW;
					dW = exvar[iex].dt_respa * exvar[iex].fcoup * exvar[iex].v;
					exvar[iex].dafed_work += 0.5*(dW + exvar[iex].dWold);
					exvar[iex].dWold = dW;
				}
				else if(parameter.useMetadyn == 2){	//metadynamics
					exvar[iex].fatom += exvar[iex].fbias;
					exvar[iex].fatom += exvar[iex].fwall;
				}
				else{
					cerr<<"\n\n\t!!ERROR!! Stupid: useMetadyn = "<<parameter.useMetadyn<<endl;
					cerr<<"\tExiting programme...\n";
					exit(1);
				}
			}
		}
		MPI_Barrier(world);
		for(iex=0;iex<nex;iex++){
			//communicate position, veolcity and force on extended variable
			MPI_Bcast(&exvar[iex].x,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].v,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].f,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fcoup,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fbias,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fwall,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fLwall,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fUwall,1,MPI_DOUBLE,0,world);
			MPI_Bcast(&exvar[iex].fatom,1,MPI_DOUBLE,0,world);
		}
		//initialize virial contribution
		for(i=0;i<6;i++){
			Fix::virial[i] = 0.0;
		}
//DK: CHECK HERE
		for(iex=0;iex<nex;iex++){
			//update the forces on all atoms
			//in UFED this is only the part from the coupling
			for(i=0;i<nlocal;i++){
				for(j=0;j<3;j++){
					//atom->f[i][j] -= exvar[iex].fcoup*exvar[iex].dQ[i][j];
					atom->f[i][j] += exvar[iex].fatom*exvar[iex].dQ[i][j];
				}
			}
			//update virial
			//in UFED this is only the part from the coupling
			//Fix::virial[0] += -exvar[iex].fcoup*exvar[iex].dvirial[0][0];				
			//Fix::virial[1] += -exvar[iex].fcoup*exvar[iex].dvirial[1][1];				
			//Fix::virial[2] += -exvar[iex].fcoup*exvar[iex].dvirial[2][2];				
			//Fix::virial[3] += -exvar[iex].fcoup*exvar[iex].dvirial[0][1];				
			//Fix::virial[4] += -exvar[iex].fcoup*exvar[iex].dvirial[0][2];				
			//Fix::virial[5] += -exvar[iex].fcoup*exvar[iex].dvirial[1][2];							

			Fix::virial[0] += exvar[iex].fatom*exvar[iex].dvirial[0][0];				
			Fix::virial[1] += exvar[iex].fatom*exvar[iex].dvirial[1][1];				
			Fix::virial[2] += exvar[iex].fatom*exvar[iex].dvirial[2][2];				
			Fix::virial[3] += exvar[iex].fatom*exvar[iex].dvirial[0][1];				
			Fix::virial[4] += exvar[iex].fatom*exvar[iex].dvirial[0][2];				
			Fix::virial[5] += exvar[iex].fatom*exvar[iex].dvirial[1][2];
		}
		//begin debug
		//of_debug_virial<<setw(7)<<step<<"  ";
		//of_debug_virial<<setprecision(4)<<fixed<<time<<"  ";
		//for(i=0;i<6;i++){
		//	of_debug_virial<<setprecision(6)<<scientific<<Fix::virial[i]<<"  ";
		//}
		//of_debug_virial<<endl;
		//end debug	

		MPI_Barrier(world);
	}
	
	//write some stuff if you are proc-0
	if(me==0){
		if(count_cv == nprint_cv){ 
			time = update->atime + (update->ntimestep-update->atimestep)*update->dt;
			print_colvar(step,time);	//print exvar/colvar values to file
			if(parameter.useNN == 1){
				print_nnout(step,time);		//print NN output to file
			}

			if(parameter.useMetadyn < 2){		//dafed and ufed
				of_ene<<setw(7)<<step<<"  ";
				of_ene<<setprecision(4)<<fixed<<time<<"  ";
				for(iex=0;iex<nex;iex++){
					print_dafed(exvar[iex]);
				}
				of_ene<<endl;
			}

			//if using ufed and metadynamics
			//print all forces
			of_biasforce<<setw(7)<<step<<"  ";
			of_biasforce<<setprecision(4)<<fixed<<time<<"  ";
			for(iex=0;iex<nex;iex++){
				print_biasforce(exvar[iex]);
			}
			of_biasforce<<endl;

			//reset counter
			count_cv = 0;
		}
//		else{	//this is only on proc-0
//			count_cv += 1;
//		}
		if(parameter.useMetadyn > 0){
			if(count_printbias == nprint_bias){
				//print ufed bias potential
				print_biaspot();
				count_printbias = 0;
			}
		}
				
	}
	//everybody writes derivatives
	if(count_deriv == nprint_deriv){
		of_deriv<<"#Time step: "<<step<<endl;
		for(i=0;i<nlocal;i++){
			of_deriv<<setw(6)<<atom->tag[i];
			for(iex=0;iex<nex;iex++){
				for(j=0;j<3;j++){
					of_deriv<<setw(18)<<setprecision(6)<<scientific<<exvar[iex].dQ[i][j];
				}
			}
			of_deriv<<endl;
		}
		count_deriv=0;
	}


	//update print counters
	if(me==0){	//this is only on proc-0
		count_cv += 1;
		if(parameter.useMetadyn > 0){
			count_printbias += 1;
		}
	}
	count_deriv += 1;


}


// ------------------------------------------------------------------------
// post force respa
// ------------------------------------------------------------------------
void FixDafedgen::post_force_respa(int vflag, int ilevel, int iloop)
{
	if (ilevel == nlevels_respa-1) post_force(vflag);
}


// ------------------------------------------------------------------------
// end of the step
// ------------------------------------------------------------------------
void FixDafedgen::end_of_step()
{
	
	int step=update->ntimestep;
	int i,j;
	int natoms=int(atom->natoms);
	int nlocal=int(atom->nlocal);
	double v_q,v_q_tot;
	double m_q,m_q_tot;
	double *mass = atom->mass;
	int *type = atom->type;
	int iex,nex;
	double m_scale;
	nex = parameter.nex;


	//determine velocity of collective variable; and the mass
	// dq/dt = \sum_i dq/dx_i * dx_i/dt = sum_i dq/dx_i * v_i
	// m_cv = 1/(\sum_i (dq/dx_i)^2/m_i)  for metal units the mass of atoms i in [g/mol]
	// => m_cv = [g/mol * Ang^2], convert to [eV * ps^2]:  1 g/mol*Ang^2 = 1.036427e-04 eV*ps^2
	//DK: CHECK HERE
	//For real units, the mass of atoms is in [g/mol], distance in [Ang] which means I can use the same conversion as metals.
    	if (force->boltz == 0.0019872067){	//for real
		m_scale = 2.390412e3;	
	}
    	else if (force->boltz == 8.617343e-5){	//for metal
		m_scale = 1.036427e-04;	
	
	}
	if(step==0){
		for(iex=0;iex<nex; iex++){
			v_q = 0.0;
			m_q = 0.0;
			for(i=0;i<nlocal;i++){
				for(j=0;j<3;j++){
					v_q += atom->v[i][j]*exvar[iex].dQ[i][j];
					m_q += exvar[iex].dQ[i][j]*exvar[iex].dQ[i][j]/mass[type[i]];
				}
				//of<<"Mass of atom "<<i<<" : "<<mass[type[i]]<<endl;
			}
			MPI_Allreduce(&m_q,&m_q_tot,1,MPI_DOUBLE,MPI_SUM,world);
			MPI_Allreduce(&v_q,&v_q_tot,1,MPI_DOUBLE,MPI_SUM,world);
			m_q_tot = 1.0/m_q_tot;

			if(me==0){
				of_colvar<<v_q_tot<<"  ";
				of_colvar<<m_q_tot*m_scale<<"  ";
			}
		}
		if(me==0){
			of_colvar<<endl;
		}
	}
	else if (count2_cv == nprint_cv){
		for(iex=0;iex<nex;iex++){
			v_q = 0.0;
			m_q = 0.0;
			for(i=0;i<nlocal;i++){
				for(j=0;j<3;j++){
					v_q += atom->v[i][j]*exvar[iex].dQ[i][j];
					m_q += exvar[iex].dQ[i][j]*exvar[iex].dQ[i][j]/mass[type[i]];
				}
			}
	
			MPI_Allreduce(&m_q,&m_q_tot,1,MPI_DOUBLE,MPI_SUM,world);
			MPI_Allreduce(&v_q,&v_q_tot,1,MPI_DOUBLE,MPI_SUM,world);	
			m_q_tot = 1.0/m_q_tot;
			if(me==0){
				of_colvar<<v_q_tot<<"  ";
				of_colvar<<m_q_tot*m_scale<<"  ";
			}
		}
		if(me==0){
			of_colvar<<endl;
		}
		
		count2_cv = 0;
	}
	count2_cv += 1;
//DK:CHECK HERE
	//copy dQ/dx derivatives to array
	
	for(i=0;i<nlocal;i++){
		for(j=0;j<3;j++){
			array[i][j] = exvar[0].dQ[i][j]; 	//hard-coded for 1 exvar
		}
		array[i][3] = exvar[0].fatom;		//total force on atom
		array[i][4] = exvar[0].fcoup;		//dafed coupling force
		array[i][5] = exvar[0].fbias;		//metadynamics bias force
		array[i][6] = exvar[0].fLwall;		//lower wall
		array[i][7] = exvar[0].fUwall;		//upper wall
		//values for restraining is set in functions get_pathcv_values_Xpoints()
		//array[i][8-10]
	}
	

	

}


// --------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------
// THE FOLLOWING ARE ALL HELPER FUNCTIONS! (INTEGRATION, ORDER PARAMETER CALCULATION, PRINTING, ETC...)
// --------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------


// ------------------------------------------------------------------------
// open colvar file and print header (in the 'constructor')
// ------------------------------------------------------------------------
void FixDafedgen::init_colvar_file(string &cvfile)
{
	//open colvar file
	of_colvar.open(cvfile.c_str(),ofstream::out | ofstream::trunc);
	
	if(parameter.choose_exvar==0){
		if(parameter.nex == 1){
			of_colvar<<"#step    time [ps]";
			of_colvar<<"     x_cv_a15 []      x_s_a15 []     v_s_a15 [1/ps]     f_s_a15 [eV]";
			of_colvar<<"     v_cv_a15 [1/ps]   m_cv_a15 [eV*ps^2]";
			of_colvar<<endl;
		}
		else if(parameter.nex == 2){
			of_colvar<<"#step    time [ps]   x_cv_bcc []      x_s_bcc []     v_s_bcc [1/ps]     f_s_bcc [eV]";
			of_colvar<<"     x_cv_a15 []      x_s_a15 []     v_s_a15 [1/ps]     f_s_a15 [eV]";
			of_colvar<<"     v_cv_bcc [1/ps]   m_cv_bcc [eV*ps^2]     v_cv_a15 [1/ps]   m_cv_a15 [eV*ps^2]";
			of_colvar<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for A15) or 2 (for BCC and A15) !");
		}
	}
	else if(parameter.choose_exvar==1){
		if(parameter.nex==1){
    			if (force->boltz == 0.0019872067){
				of_colvar<<"#step    time [ps]";
				of_colvar<<"     x_cv_f(Q) []       x_s_f(Q) []      v_s_f(Q) [1/fs]      f_s_f(Q) [kCal/mol]";
				of_colvar<<"     v_cv_f(Q) [1/fs]    m_cv_f(Q) [kCal/mol*fs^2]";
				of_colvar<<endl;
			}
    			else if (force->boltz == 8.617343e-5){
				of_colvar<<"#step    time [ps]";
				of_colvar<<"     x_cv_f(Q) []       x_s_f(Q) []      v_s_f(Q) [1/ps]      f_s_f(Q) [eV]";
				of_colvar<<"     v_cv_f(Q) [1/ps]    m_cv_f(Q) [eV*ps^2]";
				of_colvar<<endl;
			}
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for f(QBCC,QA15)) !");
		}
	}
	else if(parameter.choose_exvar==2){
		if(parameter.nex == 1){
			of_colvar<<"#step    time [ps]";
			of_colvar<<"     x_cv_Q6 []       x_s_Q6 []      v_s_Q6 [1/ps]      f_s_Q6 [eV]";
			of_colvar<<"     v_cv_Q6 [1/ps]    m_cv_Q6 [eV*ps^2]";
			of_colvar<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for Q6) !");
		}
	}
	else{
		error->all(FLERR,"Variable choose_exvar can only be 0, 1 or 2!");
	}
}

// ------------------------------------------------------------------------
// open nnout file and print header (in the 'constructor')
// ------------------------------------------------------------------------
void FixDafedgen::init_nnout_file()
{
	of_nnout.open("nnout.dat",ofstream::out | ofstream::trunc);
	
	of_nnout<<"#step    time[ps]      ureaI     ureaIV    liq    ureaA     ureaB     ureaC"<<endl;
}

// ------------------------------------------------------------------------
// open dafed_ene file and print header (in the 'constructor')
// ------------------------------------------------------------------------
void FixDafedgen::init_dafedene_file()
{
	of_ene.open("dafed_ene.dat",ofstream::out | ofstream::trunc);
	
	if(parameter.choose_exvar==0){
		if(parameter.nex == 1){
			of_ene<<"#step    time[ps]"; 
			of_ene<<"   x_s_a15 []   T_s_a15 [K]   Etot_s_a15 [eV] Wtrans_a15 [eV]"; 
			of_ene<<endl;
		}
		else if(parameter.nex == 2){
			of_ene<<"#step    time[ps]   x_s_bcc []   T_s_bcc [K]   Etot_s_bcc [eV] Wtrans_bcc [eV]"; 
			of_ene<<"   x_s_a15 []   T_s_a15 [K]   Etot_s_a15 [eV] Wtrans_a15 [eV]"; 
			of_ene<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for A15) or 2 (for BCC and A15) !");
		}
	}
	else if(parameter.choose_exvar==1){
		if(parameter.nex == 1){
			of_ene<<"#step    time[ps]"; 
			of_ene<<"   x_s_f(Q) []    T_s_f(Q) [K]    Etot_s_f(Q) [eV]  Wtrans_f(Q) [eV]"; 
			of_ene<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for f(QBCC,QA15)) !");
		}
	}
	else if(parameter.choose_exvar==2){
		if(parameter.nex == 1){
			of_ene<<"#step    time[ps]"; 
			of_ene<<"   x_s_Q6 []    T_s_Q6 [K]    Etot_s_Q6 [eV]  Wtrans_Q6 [eV]"; 
			of_ene<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for Q6) !");
		}
	}
	else{
		error->all(FLERR,"Variable choose_exvar can only be 0, 1  or 2!");
	}
}

// ------------------------------------------------------------------------
// open files to record forces and bias potential on extended variable in UFED (in the 'constructor')
// ------------------------------------------------------------------------
void FixDafedgen::init_bias_files()
{
	of_biasforce.open("all_forces.dat",ofstream::out | ofstream::trunc);
	of_biaspot.open("ufed_potential.dat",ofstream::out | ofstream::trunc);
			
	if(parameter.choose_exvar==0){
		if(parameter.nex == 1){
			of_biasforce<<"#step    time[ps]";
			of_biasforce<<"   f_s_a15 [eV]   fatom_s_a15 [eV]   fcoup_s_a15 [eV]   fbias_s_a15 [eV]   fwall_s_a15 [eV]   fLwall_s_a15 [eV]   fUwall_s_a15 [eV]";
			of_biasforce<<endl;

			of_biaspot<<"#s_a15   Vbias(s)";
			of_biaspot<<endl;
		}
		else if(parameter.nex==2){
			of_biasforce<<"#step    time[ps]";
			of_biasforce<<"   f_s_bcc [eV]   fatom_s_bcc [eV]   fcoup_s_bcc [eV]   fbias_s_bcc [eV]   fwall_s_bcc [eV]   fLwall_s_bcc [eV]   fUwall_s_bcc [eV]";
			of_biasforce<<"   f_s_a15 [eV]   fatom_s_a15 [eV]   fcoup_s_a15 [eV]   fbias_s_a15 [eV]   fwall_s_a15 [eV]   fLwall_s_a15 [eV]   fUwall_s_a15 [eV]";
			of_biasforce<<endl;
					
			of_biaspot<<"#s_bcc   s_a15   Vbias(s)";
			of_biaspot<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for A15) or 2 (for BCC and A15) !");
		}
	}
	else if(parameter.choose_exvar==1){
		if(parameter.nex == 1){
			of_biasforce<<"#step    time[ps]";
			if(parameter.use_restraint==1){
				of_biasforce<<"   f_s_f(Q) [eV]   fatom_s_f(Q) [eV]   fcoup_s_f(Q) [eV]   fbias_s_f(Q) [eV]   fwall_s_f(Q) [eV]   fLwall_s_f(Q) [eV]   fUwall_s_f(Q) [eV]   frestraint_z [eV]      vrestraint_z [eV]";
			}
			else{
				of_biasforce<<"   f_s_f(Q) [eV]   fatom_s_f(Q) [eV]   fcoup_s_f(Q) [eV]   fbias_s_f(Q) [eV]   fwall_s_f(Q) [eV]   fLwall_s_f(Q) [eV]   fUwall_s_f(Q) [eV]";
			}

			of_biasforce<<endl;

			of_biaspot<<"#s_f(Q)   Vbias(s)";
			of_biaspot<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for f(QBCC,QA15)) !");
		}
	}
	else if(parameter.choose_exvar==2){
		if(parameter.nex == 1){
			of_biasforce<<"#step    time[ps]";
			of_biasforce<<"   f_s_Q6 [eV]   fatom_s_Q6 [eV]   fcoup_s_Q6 [eV]   fbias_s_Q6 [eV]   fwall_s_Q6 [eV]   fLwall_s_Q6 [eV]   fUwall_s_Q6 [eV]";
			of_biasforce<<endl;

			of_biaspot<<"#s_Q6   Vbias(s)";
			of_biaspot<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for Q6) !");
		}
	}
	else{
		error->all(FLERR,"Variable choose_exvar can only be 0, 1 or 2!");
	}
}

// ------------------------------------------------------------------------
// open derivatives files and print headers, one for each proc (in the 'constructor')
// ------------------------------------------------------------------------
void FixDafedgen::init_derivatives_files(string &dcvfile)
{
	ostringstream ss;
	ss<<me;
	dcvfile = dcvfile+"-"+ss.str();
	
	//OPEN ONE FILE FOR EACH PROC
	of_deriv.open(dcvfile.c_str(),ofstream::out | ofstream::trunc);
	if(parameter.choose_exvar==0){
		if(parameter.nex==1){
			of_deriv<<"#atom_id";   
			of_deriv<<"      dq_a15/dx      dq_a15/dy      dq_a15/dz"<<endl;   
		}
		else if (parameter.nex==2){
			of_deriv<<"#atom_id      dq_bcc/dx      dq_bcc/dy      dq_bcc/dz";
			of_deriv<<"      dq_a15/dx      dq_a15/dy      dq_a15/dz"<<endl;
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for A15) or 2 (for BCC and A15) !");
		}
	}
	else if(parameter.choose_exvar==1){
		if(parameter.nex == 1){
			of_deriv<<"#atom_id";   
			of_deriv<<"      dq_f(Q)/dx       dq_f(Q)/dy       dq_f(Q)/dz"<<endl;   
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for f(QBCC,QA15)) !");
		}
	}
	else if(parameter.choose_exvar==2){
		if(parameter.nex == 1){
			of_deriv<<"#atom_id";   
			of_deriv<<"      dq_Q6/dx       dq_Q6/dy       dq_Q6/dz"<<endl;   
		}
		else{
			error->all(FLERR,"Number of collective variables can only be 1 (for Q6) !");
		}
	}
	else{
		error->all(FLERR,"Variable choose_exvar can only be 0, 1 or 2!");
	}

}


// ------------------------------------------------------------------------
// GGMT integrator for the extended variable
// ------------------------------------------------------------------------
void FixDafedgen::integrate_ggmt(DAFED::EXvariables &exvar)
{
	//int step=update->ntimestep;
	DAFED::CGgmt *th = &(exvar.ggmt);
	int iii, jjj;
	double aa=0.0, bb=0.0;
	double kt;
	double G1,G2;
	double dt2,dt4,dt8;

	// mv^2
	kt = exvar.m*exvar.v*exvar.v;

	//Loop over the RESPA integration steps
	for(iii=1; iii<=th->n_respa_ggmt; iii++){
		//Loop over the three terms of the Suzuki-Yoshida decomposition:
		for(jjj=1; jjj<=3; jjj++){
			
			dt2 = 0.5*th->dt[jjj];		// half time step
			dt4 = 0.25*th->dt[jjj];		// 1/4 time step
			dt8 = 0.125*th->dt[jjj];	// 1/8 time step

			G1 = (kt - th->kT0)/(th->Q1);
			G2 = (kt*kt/3.0 - th->kT0*th->kT0)/(th->Q2);
			th->v1 += dt4*G1;
			th->v2 += dt4*G2;

			aa = exp(-dt8*(th->v1 + th->kT0*th->v2));
			exvar.v *= aa;
			kt = exvar.m*exvar.v*exvar.v;

			bb = kt*th->v2/3.0;
			exvar.v *= sqrt(1.0/(1.0 + dt2 * bb));

			exvar.v *= aa;
			kt = exvar.m*exvar.v*exvar.v;

			th->eta1 += dt2*th->v1;
			th->eta2 += dt2*th->v2*(th->kT0 + kt);

			aa = exp(-dt8*(th->v1 + th->kT0*th->v2));
			exvar.v *= aa;
			kt = exvar.m*exvar.v*exvar.v;

			bb = kt*th->v2/3.0;
			exvar.v *= sqrt(1.0/(1.0 + dt2 * bb));

			exvar.v *= aa;
			kt = exvar.m*exvar.v*exvar.v;

			G1 = (kt - th->kT0)/(th->Q1);
			G2 = (kt*kt/3.0 - th->kT0*th->kT0)/(th->Q2);
			th->v1 += dt4*G1;
			th->v2 += dt4*G2;
		}
	}
}


// ------------------------------------------------------------------------
// print colvar/exvar to file
// ------------------------------------------------------------------------
void FixDafedgen::print_colvar(int step, double time)
{
	int iex;
	int nex=parameter.nex;

	of_colvar<<setw(7)<<step<<"  ";
	of_colvar<<setprecision(4)<<fixed<<time<<"  ";
	for(iex=0;iex<nex;iex++){
		of_colvar<<setprecision(6)<<scientific<<exvar[iex].Q<<"  ";
		of_colvar<<exvar[iex].x<<"  "<<exvar[iex].v<<"  "<<exvar[iex].f<<"  ";
	}
}

// ------------------------------------------------------------------------
// print NN out put to file
// ------------------------------------------------------------------------
void FixDafedgen::print_nnout(int step, double time)
{
	int iex;
	of_nnout<<setw(7)<<step<<"  ";
	of_nnout<<setprecision(4)<<fixed<<time<<"  ";
//DK: CHECK HERE is iex really 6? <-- YES. You have 6 outputs and you defined Qall as 6 in fix_dafed.h file
	for(iex=0;iex<6;iex++){
		of_nnout<<setprecision(6)<<scientific<<Qall[iex]<<"  ";
	}
	of_nnout<<endl;
}

// ------------------------------------------------------------------------
// print some energy stuff from the dAFED calculations
// ------------------------------------------------------------------------
void FixDafedgen::print_dafed(DAFED::EXvariables &exvar)
{
	double potential;
	double energy;
	double temperature;
	double transferred_work;
	double delta;
	double energy_ggmt;

	delta = exvar.Q - exvar.x;

	potential = 0.5*exvar.kappa*delta*delta;
	
	energy_ggmt = 0.5*(exvar.ggmt.Q1*exvar.ggmt.v1*exvar.ggmt.v1) +
		0.5*(exvar.ggmt.Q2*exvar.ggmt.v2*exvar.ggmt.v2) +
		exvar.ggmt.kT0*(exvar.ggmt.eta1 + exvar.ggmt.eta2);
	
	energy = 0.5*exvar.m*exvar.v*exvar.v + potential + energy_ggmt;

	temperature = exvar.v*exvar.v*exvar.m/force->boltz;

	// What we want is not the accumulated work which pertains to the system H_phys + V
	// But we want the work W transferred to the physical system, which pertains to H_phys
	// Where H_phys is the energy of the physical system + thermostat
	// and V is the potential 0.5*kappa*(q-s)^2
	// For a discussion of this (in the context of steered MD), see
	// Schurr and Fujimoto, J. Phys. Chem. B 107. 14007 (2003)
	// This way, both H-W and H_dafed+W should be conserved
	transferred_work = exvar.dafed_work - potential;

	of_ene<<setprecision(6)<<scientific<<exvar.x<<"  ";
	of_ene<<setprecision(6)<<scientific<<temperature<<"  ";
	of_ene<<setprecision(6)<<scientific<<energy<<"  ";
	of_ene<<setprecision(6)<<scientific<<transferred_work<<"  ";

}

//------------------------------------------------------------------------
// print force due to ufed bias potential
//------------------------------------------------------------------------
void FixDafedgen::print_biasforce(DAFED::EXvariables &exvar){
	
	of_biasforce<<setprecision(6)<<scientific<<exvar.f<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.fatom<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.fcoup<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.fbias<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.fwall<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.fLwall<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.fUwall<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.frestrain<<"  ";
	of_biasforce<<setprecision(6)<<scientific<<exvar.vrestrain<<"  ";
}

//------------------------------------------------------------------------
// print accumulated ufed bias potential
//------------------------------------------------------------------------
void FixDafedgen::print_biaspot(){
	
	int i,iex;
	int nex=parameter.nex;

	for(i=0;i<Nlattot;i++){
		for(iex=0;iex<nex;iex++){
			of_biaspot<<setprecision(6)<<scientific<<bhist[i].x[iex]<<"  ";
		}
		of_biaspot<<setprecision(6)<<scientific<<bhist[i].pot;
		of_biaspot<<endl;
	}
	of_biaspot<<endl;


}

//------------------------------------------------------------------------
//  calculate the order parameter in two steps, since we need information
//  of neighbours to calculate the derivatives of the Steinhardt parameters
//------------------------------------------------------------------------
/*void FixDafed::get_exvar_values()
{
	int i,j;
	int L,mi;
	double box[3];				//size of simulation box
	double **x = atom->x;			//atom vector, contains only local atoms, owned and ghost
	tagint *tag = atom->tag;
 
	
	// only orthorhombic box works here!
	box[0] = domain->h[0];
	box[1] = domain->h[1];
	box[2] = domain->h[2];

	//initialize data structure that hold info on lammps neighbour list
	lmpneigh.inum = list->inum;
	lmpneigh.ilist = list->ilist;
	lmpneigh.numneigh = list->numneigh;
	lmpneigh.firstneigh = list->firstneigh;
	lmpneigh.nall = atom->nlocal + atom->nghost;
	lmpneigh.me = me;
	lmpneigh.ntot = atom->natoms;

	//copy current positions into data structure
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			molecules[i].pos[j] = x[i][j];
		}
	}
	
	//make small neighbour lists for calculating collective variable 
	DAFED::get_symmetry_NN_3_lmpneighbour(molecules,lmpneigh,box,parameter,sf2Rs,sf2eta,sf3kappa,NNvar);


	//from fcc/orient: copy stuff into nbr structure and communicate
	for(i=0;i<atom->nlocal;i++){
		for(j=0;j<molecules[i].n_neighbors1;j++){
			nbr[i].id[j] = molecules[i].neighbors1[j];
		}
	}

	// communicate to acquire nbr and molecules data for ghost atoms
	comm->forward_comm_fix(this);

	//call second part after communication
	//this is specifically for nex=2 with 0=BCC and 1=A15
	//and                  for nex=1 with 0=A15
	DAFED::get_Q_3_lmpneighbour(molecules,lmpneigh,exvar,parameter,sf2Rs,sf2eta,sf3kappa,NNvar,Qlocal);
	MPI_Barrier(world);
 
	//sum contributions from all procs
	for(i=0;i<parameter.nex;i++){
		MPI_Allreduce(&exvar[i].Qlocal,&exvar[i].Q,1,MPI_DOUBLE,MPI_SUM,world);
	}
	MPI_Allreduce(&Qlocal,&Qall,5,MPI_DOUBLE,MPI_SUM,world);
}
*/
//------------------------------------------------------------------------
//  calculate the order parameter in two steps, since we need information
//  of neighbours to calculate the derivatives of the Steinhardt parameters
//  get path collective variable as a function of QBCC and QA15
//------------------------------------------------------------------------
/*void FixDafed::get_pathcv_values()
{
	int i,j;
	int iex,nex;
	int L,mi;
	double box[3];				//size of simulation box
	double **x = atom->x;			//atom vector, contains only local atoms, owned and ghost
	tagint *tag = atom->tag;
	DAFED::EXvariables exvar_tmp[2];

	//initialize some stuff for temporary exvar
	nex = 2;
	for(iex=0;iex<nex;iex++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				exvar_tmp[iex].dvirial[i][j] = 0.0;
			}
		}
	}
	//allocate memory for extended variable derivatives
	for(iex=0;iex<nex;iex++){
		//allocate derivatives
		exvar_tmp[iex].dQ = new double* [nmax];
		for(int j=0;j<nmax;j++){
			exvar_tmp[iex].dQ[j] = new double [3];
		}
	}

	
	// only orthorhombic box works here!
	box[0] = domain->h[0];
	box[1] = domain->h[1];
	box[2] = domain->h[2];

	//initialize data structure that hold info on lammps neighbour list
	lmpneigh.inum = list->inum;
	lmpneigh.ilist = list->ilist;
	lmpneigh.numneigh = list->numneigh;
	lmpneigh.firstneigh = list->firstneigh;
	lmpneigh.nall = atom->nlocal + atom->nghost;
	lmpneigh.me = me;
	lmpneigh.ntot = atom->natoms;

	//copy current positions into data structure
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			molecules[i].pos[j] = x[i][j];
		}
	}
	
	//make small neighbour lists for calculating collective variable 
	DAFED::get_symmetry_NN_3_lmpneighbour(molecules,lmpneigh,box,parameter,sf2Rs,sf2eta,sf3kappa,NNvar);


	//from fcc/orient: copy stuff into nbr structure and communicate
	for(i=0;i<atom->nlocal;i++){
		for(j=0;j<molecules[i].n_neighbors1;j++){
			nbr[i].id[j] = molecules[i].neighbors1[j];
		}
	}

	// communicate to acquire nbr and molecules data for ghost atoms
	comm->forward_comm_fix(this);

	//call second part after communication
	//this is specifically for nex=2 with 0=BCC and 1=A15
	//and                  for nex=1 with 0=A15
	//set parameter.nex = 2 since we want both QBCC and QA15
	parameter.nex = 2;
	DAFED::get_Q_3_lmpneighbour(molecules,lmpneigh,exvar_tmp,parameter,sf2Rs,sf2eta,sf3kappa,NNvar,Qlocal);
	//reset parameter.nex to 1 since only a singel path collective variable
	parameter.nex = 1;
	MPI_Barrier(world);

	//sum contributions from all procs
	for(i=0;i<nex;i++){
		MPI_Allreduce(&exvar_tmp[i].Qlocal,&exvar_tmp[i].Q,1,MPI_DOUBLE,MPI_SUM,world);
	}
	MPI_Allreduce(&Qlocal,&Qall,5,MPI_DOUBLE,MPI_SUM,world);

	
	//here we have now the values of Q-BCC, Q-A15 as well as all dQ-BCC/dx and dQ-A15/dx values
	double QBCC_1, QBCC_2;		//values in state 1 and 2
	double QA15_1, QA15_2;
	double D1, D2;		//distance D1 and D2 in terms of QBCC and QA15 values;
	double expD1, expD2;
	double lambda;
	double pathcv;
	double dpcv_dBCC, dpcv_dA15;

	QBCC_1 = 0.2;		//state 1: 0.2-BCC and 0.6-A15
	QBCC_2 = 0.6;		//state 2: 0.6-BCC and 0.2-A15
	QA15_1 = 0.6;
	QA15_2 = 0.2;

	lambda = 10.0;

	//exvar_tmp[0] = BCC, exvar_tmp[1] = A15
	D1 = sqrt((QBCC_1-exvar_tmp[0].Q)*(QBCC_1-exvar_tmp[0].Q) + (QA15_1-exvar_tmp[1].Q)*(QA15_1-exvar_tmp[1].Q));
	D2 = sqrt((QBCC_2-exvar_tmp[0].Q)*(QBCC_2-exvar_tmp[0].Q) + (QA15_2-exvar_tmp[1].Q)*(QA15_2-exvar_tmp[1].Q));

	expD1 = exp(-lambda*D1);
	expD2 = exp(-lambda*D2);

	pathcv = (expD1 + 2.0*expD2)/(expD1 + expD2); 

	//assign value to cv used in simulation
	exvar[0].Q = pathcv;

	//now get the derivative as dS/dx = dS/dBCC*dBCC/dx + dS/dA15*dA15/dx
	dpcv_dBCC = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QBCC_1-exvar_tmp[0].Q)/D1 + expD2*lambda*(QBCC_2-exvar_tmp[0].Q)/D2))
			/((expD1 + expD2)*(expD1 + expD2)))
		+ (expD1*lambda*(QBCC_1-exvar_tmp[0].Q)/D1 + 2.0*expD2*lambda*(QBCC_2-exvar_tmp[0].Q)/D2)/(expD1 + expD2);

	dpcv_dA15 = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QA15_1-exvar_tmp[1].Q)/D1 + expD2*lambda*(QA15_2-exvar_tmp[1].Q)/D2))
			/((expD1 + expD2)*(expD1 + expD2)))
		+ (expD1*lambda*(QA15_1-exvar_tmp[1].Q)/D1 + 2.0*expD2*lambda*(QA15_2-exvar_tmp[1].Q)/D2)/(expD1 + expD2);

	for(i=0;i<atom->nlocal;i++){
		for(j=0;j<3;j++){
			exvar[0].dQ[i][j] = dpcv_dBCC * exvar_tmp[0].dQ[i][j] + dpcv_dA15 * exvar_tmp[1].dQ[i][j];
		}
	}

	//and also the virial (this is the local sum for each proc, also in lammps these are determined for each proc)
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			exvar[0].dvirial[i][j] = (dpcv_dBCC * exvar_tmp[0].dvirial[i][j]) + (dpcv_dA15 * exvar_tmp[1].dvirial[i][j]);
		}
	}

	//check if additional restraint is applied
	//
	//  !! THIS IS THE OLD RESTRAIN ON THE SUM of Qbcc + QA15 !!
	//
	if(parameter.use_restraint == 1){
		double Qsum;
		double diffrest, diffrest_prod;

		Qsum = exvar_tmp[0].Q + exvar_tmp[1].Q;
		diffrest = Qsum - parameter.restraint_value;
		diffrest_prod = 1.0;
		for(int iwall=0; iwall < (parameter.restraint_n-1); iwall++){
			diffrest_prod *= diffrest;
		}
		exvar_tmp[0].frestrain = parameter.restraint_pref * diffrest_prod;	//-dVres/dq
		exvar_tmp[1].frestrain = parameter.restraint_pref * diffrest_prod;
		exvar[0].frestrain = parameter.restraint_pref * diffrest_prod;
		
		for(iex=0;iex<nex;iex++){
			//add forces to atom
			for(i=0;i<atom->nlocal;i++){
				for(j=0;j<3;j++){
					atom->f[i][j] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dQ[i][j];
				}
			}
			//update virial
			Fix::virial[0] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dvirial[0][0];				
			Fix::virial[1] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dvirial[1][1];				
			Fix::virial[2] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dvirial[2][2];				
			Fix::virial[3] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dvirial[0][1];				
			Fix::virial[4] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dvirial[0][2];				
			Fix::virial[5] += exvar_tmp[iex].frestrain*exvar_tmp[iex].dvirial[1][2];
		}
	}

	
	//free memory of derivatives of extended variables
	for(iex=0;iex<nex;iex++){
		for(j=0;j<nmax;j++){
			delete [] exvar_tmp[Uliqiex].dQ[j];
		}
		delete [] exvar_tmp[iex].dQ;
	}


}*/

void FixDafedgen::get_pathcv_values_gen()
{
	int i,j;
	int iex,nex;
	int L,mi;
	double box[6];				//size of simulation box
	double **x = atom->x;			//atom vector, contains only local atoms, owned and ghost
	tagint *tag = atom->tag;
	tagint *moltag = atom->molecule;
	int *type = atom->type;
	DAFED::EXvariables exvar_tmp[2];
	int natoms=int(atom->natoms);
	nex = 2;

	int jnum;
	int *jlist;
	for(iex=0;iex<nex;iex++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				exvar_tmp[iex].dvirial[i][j] = 0.0;
			}
		}
	}
	
	cout << "lmpneigh.nall  = " << lmpneigh.nall << endl;
	cout << "nmax = " << nmax << endl;

	//allocate memory for extended variable derivatives
	for(iex=0;iex<nex;iex++){
		//allocate derivatives
		exvar_tmp[iex].dQ = new double* [nmax];
		for(int j=0;j<nmax;j++){
			exvar_tmp[iex].dQ[j] = new double [3];
		}
	}
	//NNout is used for the trajectory data output. 
	for (i = 0;i<nmax;i++)        //loop over a maximum number of atoms
	{
		for (int outi = 0; outi < parameter.nnout; outi++)    //loop over number of NNoutput
		{                                       //set up values to be zero for NNout 
			molecules[i].NNout[outi] = 0;
			NNout[i][outi] = 0;
		}//end of NNout loop
	}//end of nmax loop	
	// only orthorhombic box works here!
	parameter.triclinic = domain->triclinic;
	box[0] = domain->h[0];
	box[1] = domain->h[1];
	box[2] = domain->h[2];
	if (parameter.triclinic)
	{
		box[3] = domain->h[3];//yz
        	box[4] = domain->h[4];//xz
       		box[5] = domain->h[5];//xy
	}
	else
	{
		box[3] = 0;
        	box[4] = 0;
       		box[5] = 0;	
	}
	//initialize data structure that hold info on lammps neighbour list
	lmpneigh.inum = list->inum;			//same as local num except using  pair_style hybrid. In this case, number of neigh
	lmpneigh.ilist = list->ilist;			//get local id for this specific ato
	lmpneigh.numneigh = list->numneigh;		//number of neighbors
	lmpneigh.firstneigh = list->firstneigh;		//list of neighbors for one specific atom
	lmpneigh.nall = atom->nlocal + atom->nghost;	//number of atoms in the processor(includes ghost)
	lmpneigh.me = me;				//proc number
	lmpneigh.ntot = atom->natoms;			//number of total atoms
	
	
	cout << "lmpneigh.nall: " << lmpneigh.nall << "  lmpneigh.inum: " << lmpneigh.inum << endl;

	for (int ii = 0; ii< lmpneigh.nall; ii++)
	{
		jlist = lmpneigh.firstneigh[ii];		//list of ti neighbors 
		jnum = lmpneigh.numneigh[ii];			//number of ti neighbors for ti
		for (int jj = 0; jj<jnum;jj++ )
		{
			int tj = jlist[jj];	
			tj &= MY_NEIGHMASK;
//			cout <<"ii = " << ii << "  neighbor# = " << jj << "  tj " << tj << endl;
		}

	}
//	exit(1);



	//copy current positions into data structure
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			molecules[i].pos[j] = x[i][j];	//get position for each atom
		}
		molecules[i].moltype = moltag[i]-1;	//get mol type for each atom
		molecules[i].eletype = type[i];		//get element type for each atom
		cout << " molecule ID: " << i << "  molecules[i].eletype = " << molecules[i].eletype  << "  molecules[i].moltype = " << molecules[i].moltype  << "  x = " << molecules[i].pos[0]<< endl;
	}
	cout << "lmpneigh.nall = " << lmpneigh.nall << endl;
	cout << "lmpneigh.ntot = " << lmpneigh.ntot << endl;	
	
//	exit(1);

	//get all the neighbor info to molecule class.
	DAFED::get_allNeighbourDistances_lmp_new(molecules,molec,lmpneigh,lmpneigh.nall,box,parameter);
	//comm->forward_comm_fix(this);//Need to communicate here to calculate COvec and NNvec
	comm->forward_comm(this);//Need to communicate here to calculate COvec and NNvec
//	MPI_Barrier(world);
//	exit(1);	
//	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
//	cout << "After neighbour " << ctime(&timenow) << endl;
	
	//get symmetry fuction and Neural Network result from here. 
	//DAFED::get_vec_symmetry_NN_3_lmpneighbour_new(molecules, molec,lmpneigh,box,parameter,Rskappa,eta,NNvar);
	DAFED::get_vec_symmetry_NN_3_lmpneighbour_new(molecules, molec,lmpneigh,box,parameter,parameter.RskappaLst,parameter.etaLst,NNvar);
	
//	exit(1);	
//	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "After get_vec_symmetry" << endl;
	
	//comm->forward_comm_fix(this);//Need to communicate here to calculate NNgrad
	comm->forward_comm(this);//Need to communicate here to calculate NNgrad
//	MPI_Barrier(world);
	cout << "After get_vec_symmetry and comm" << endl;
	
	
	//fwd all derivative info from here.
	parameter.nex = 2;
	//DAFED::get_Q_vec_lmpneigh(molecules,molec, lmpneigh,exvar_tmp,parameter,Rskappa,eta,NNvar,Qlocal);
	DAFED::get_Q_vec_lmpneigh(molecules,molec, lmpneigh,exvar_tmp,parameter,parameter.RskappaLst,parameter.etaLst,NNvar,Qlocal);
	//exit(1);	
//	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "After get Q " <<  endl;
//	exit(1);	
	parameter.nex = 1;
	MPI_Barrier(world);

	//sum contributions from all procs. get Q.
	for(i=0;i<nex;i++){
		MPI_Allreduce(&exvar_tmp[i].Qlocal,&exvar_tmp[i].Q,1,MPI_DOUBLE,MPI_SUM,world);
	}
	MPI_Allreduce(&Qlocal,&Qall,parameter.nnout,MPI_DOUBLE,MPI_SUM,world);
	
	
	double QU1[10],QU4[10];		//CHECK THAT THIS IS CONSISTENT WITH NO. OF POINTS
	double DISTQ[10];
	double expD[10];
	double lambda;			//setted up as 200
	double sumN,sumD;		//sum of numerator part in S(Y) and sum of denominator part in S(Y).
	double pathcv;			//S(Y)
	double dpcv_dU1, dpcv_dU4;
	int pathpts;			//No. of points along pathcv
//kokodemo dame
	//set how many points along the pathcv
	pathpts = 10;

	//define images along the path, state 0 is 0.40-U1 and 0.60-U4, last state is 0.80-U1, 0.20-U4 with n images
	for(i=0;i<pathpts;i++){
		//QU1[i] = 0.40 + i*0.04;
		//QU4[i] = 0.55 - i*0.04;
		//looks like U1 is less stable, so try to increase U4 instead.
		//QU1[i] = 0.96 - i*0.07;
		//QU4[i] = 0.14 + i*0.07;
		//QU1[i] = 0.8 - i*0.06;
		//QU4[i] = 0.2 + i*0.06;
		//QU1[i] = 0.40 - i*0.03;
		//QU4[i] = 0.55 + i*0.03;
		QU1[i] = 0.9 - i*0.07;
		QU4[i] = 0.1 + i*0.07;
	}

	lambda = 200.0;

	//exvar_tmp[0] = BCC, exvar_tmp[1] = A15
	//D1 = sqrt((QU1_1-exvar_tmp[0].Q)*(QU1_1-exvar_tmp[0].Q) + (QU4_1-exvar_tmp[1].Q)*(QU4_1-exvar_tmp[1].Q));
	//D2 = sqrt((QU1_2-exvar_tmp[0].Q)*(QU1_2-exvar_tmp[0].Q) + (QU4_2-exvar_tmp[1].Q)*(QU4_2-exvar_tmp[1].Q));
	for(i=0;i<pathpts;i++){
		//DISTQ[i] = sqrt((QU1[i]-exvar_tmp[0].Q)*(QU1[i]-exvar_tmp[0].Q) + (QU4[i]-exvar_tmp[1].Q)*(QU4[i]-exvar_tmp[1].Q));
		DISTQ[i] = (QU1[i]-exvar_tmp[0].Q)*(QU1[i]-exvar_tmp[0].Q) + (QU4[i]-exvar_tmp[1].Q)*(QU4[i]-exvar_tmp[1].Q);	//it is more likely diffQ^2. diffQ=sqrt(diffQU1**2 + diffQU4**2)
		expD[i] = exp(-lambda*DISTQ[i]);		//denominator part(except sum)
	}

	sumN = 0.0;
	sumD = 0.0;
	for(i=0;i<pathpts;i++){
		sumN += i*expD[i];
		sumD += expD[i];
	}
	pathcv = 1.0/(pathpts-1.0) * sumN/sumD;	//S(Y)

	//assign value to cv used in simulation
//	exvar[1].Q = pathcv;		//This is not in tmp!!!
	exvar[0].Q = pathcv;		//This is not in tmp!!!
	//cout << "pathcv = " << pathcv << endl; 

	//now get the derivative as dS/dx = dS/dU1*dU1/dx + dS/dU4*dU4/dx
	//dpcv_dU1 = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QU1_1-exvar_tmp[0].Q)/D1 + expD2*lambda*(QU1_2-exvar_tmp[0].Q)/D2))
	//		/((expD1 + expD2)*(expD1 + expD2)))
	//	+ (expD1*lambda*(QU1_1-exvar_tmp[0].Q)/D1 + 2.0*expD2*lambda*(QU1_2-exvar_tmp[0].Q)/D2)/(expD1 + expD2);

	//dpcv_dU4 = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QU4_1-exvar_tmp[1].Q)/D1 + expD2*lambda*(QU4_2-exvar_tmp[1].Q)/D2))
	//		/((expD1 + expD2)*(expD1 + expD2)))
	//	+ (expD1*lambda*(QU4_1-exvar_tmp[1].Q)/D1 + 2.0*expD2*lambda*(QU4_2-exvar_tmp[1].Q)/D2)/(expD1 + expD2);
	double dsumN_dU1,dsumD_dU1;
	double dsumN_dU4,dsumD_dU4;
	dsumN_dU1 = 0.0;
	dsumD_dU1 = 0.0;
	dsumN_dU4 = 0.0;
	dsumD_dU4 = 0.0;
	for(i=0;i<pathpts;i++){
		dsumN_dU1 += i*lambda*2.0*(QU1[i]-exvar_tmp[0].Q)*expD[i];
		dsumD_dU1 += lambda*2.0*(QU1[i]-exvar_tmp[0].Q)*expD[i];

		dsumN_dU4 += i*lambda*2.0*(QU4[i]-exvar_tmp[1].Q)*expD[i];
		dsumD_dU4 += lambda*2.0*(QU4[i]-exvar_tmp[1].Q)*expD[i];
	}

	dpcv_dU1 = 1.0/(pathpts-1.0)*(dsumN_dU1/sumD - sumN*dsumD_dU1/(sumD*sumD));
	dpcv_dU4 = 1.0/(pathpts-1.0)*(dsumN_dU4/sumD - sumN*dsumD_dU4/(sumD*sumD));

	//for(i=0;i<atom->nlocal;i++){
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			exvar[0].dQ[i][j] = dpcv_dU1 * exvar_tmp[0].dQ[i][j] + dpcv_dU4 * exvar_tmp[1].dQ[i][j];
		}
	}
	//DK:added for NNout.  
	for (int mi = 0; mi < lmpneigh.lmpmol; mi++)  //loop over lmpmol which is a number of the molecule in the proc
	{
		for (int atmi = 0; atmi <parameter.natm; atmi++)     //loop over a number of the atoms in each molecule 
        	{
                	for (int outi = 0; outi < parameter.nnout; outi++){   //loop over number of NNout
                        	molecules[molec[mi].atom_ID[atmi]].NNout[outi] = molecules[molec[mi].atom_ID[0]].NNout[outi];
//				if (mi == 0)
//				{
//					cout << "molecules[molec[mi].atom_ID[atmi]].NNout[outi] = " << molecules[molec[mi].atom_ID[atmi]].NNout[outi] << endl;
//				}
	                }//end of NNout
		}//end of atmi loop
	}//end of lmpmol loop
	//comm->reverse_comm_fix(this);	//get all the ghost info to real atoms. 
	comm->reverse_comm(this);	//get all the ghost info to real atoms. 
	//DK:Trnasfering all the NNout info to the atom_array. 
	for (int ti = 0; ti < atom->nlocal; ti++)   //loop over all owned atoms
	{
		for (i = 0; i < parameter.nnout; i++)         //loop over number of NNout
		{
			NNout[ti][i] =molecules[ti].NNout[i];
		}//end of NNout loop
	}//end of owned atoms loop




	//and also the virial (this is the local sum for each proc, also in lammps these are determined for each proc)
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			exvar[0].dvirial[i][j] = (dpcv_dU1 * exvar_tmp[0].dvirial[i][j]) + (dpcv_dU4 * exvar_tmp[1].dvirial[i][j]);
			//cout << exvar[0].dvirial[i][j] << endl;
		}
	}

	//check if additional restraint is applied (here this is the transition tube)
	cout <<"end of _gen code " << endl;
//**********************************Most likely there is a bug in this part. Be careful****************************	
	if(parameter.use_restraint == 1){
		//cout << "parameter.use_restraint == 1 " << endl;
		lambda = 1000.0;		//larger lambda value
		double dist_pathcv;	//"distance" from pathcv: -1/lambda*Log(sumD)*30
		double diffrest, diffrest_prod;
		double ddist_pathcv_dU1, ddist_pathcv_dU4;

		sumD = 0.0;
		for(i=0;i<pathpts;i++){
			expD[i] = exp(-lambda*DISTQ[i]);
			sumD += expD[i];
		}

		dist_pathcv = -30.0/lambda*log(sumD);


		diffrest = dist_pathcv - parameter.restraint_value;
		diffrest_prod = 1.0;
		for(int iwall=0; iwall < (parameter.restraint_n-1); iwall++){
			diffrest_prod *= diffrest;
		}
		//exvar_tmp[0].frestrain = parameter.restraint_pref * diffrest_prod;	//-dVres/dq
		//exvar_tmp[1].frestrain = parameter.restraint_pref * diffrest_prod;
		exvar[0].frestrain = parameter.restraint_pref * diffrest_prod;		//-dVres(s)/ds  where s=dist_pathcv
		exvar[0].vrestrain = -exvar[0].frestrain * diffrest / parameter.restraint_n;	// Vres(s) 

		//get derivative of dist_pathcv wrt U1 and U4
		dsumD_dU1 = 0.0;
		dsumD_dU4 = 0.0;
		for(i=0;i<pathpts;i++){
			dsumD_dU1 += lambda*2.0*(QU1[i]-exvar_tmp[0].Q)*expD[i];
			dsumD_dU4 += lambda*2.0*(QU4[i]-exvar_tmp[1].Q)*expD[i];
		}

		ddist_pathcv_dU1 = -30.0/(lambda*sumD)*dsumD_dU1;
		ddist_pathcv_dU4 = -30.0/(lambda*sumD)*dsumD_dU4;

		//start: this is just for testing the derivative
		//for(i=0;i<atom->nlocal;i++){
		//	for(j=0;j<3;j++){
		//		exvar[0].dQ[i][j] = -exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dA15*exvar_tmp[1].dQ[i][j]);
		//	}
		//}
		
		//exvar[0].Q = -parameter.restraint_pref/parameter.restraint_n*diffrest_prod*diffrest;

		//end: only for testing derivatives
		
		//add forces to atom
//		cout << "Before adding force to the atom " << endl;
//		DK: CHECK HERE AGAIN. Add forces to atoms crash the code. 
		for(i=0;i<atom->nlocal;i++){
			for(j=0;j<3;j++){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//DK:Most Likely THIS IS THE PROBLEM When use rstraint is 1
				atom->f[i][j] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dU4*exvar_tmp[1].dQ[i][j]);
				//to print out forces due to restrain
				//0-7 are the other forces!
				array[i][8+j] = exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dU4*exvar_tmp[1].dQ[i][j]);
			}
		}
//		cout << "After adding force " << endl;
//DK:commented out for now		//update virial
		Fix::virial[0] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[0][0] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[0][0]);		
		Fix::virial[1] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[1][1] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[1][1]);		
		Fix::virial[2] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[2][2] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[2][2]);		
		Fix::virial[3] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[0][1] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[0][1]);		
		Fix::virial[4] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[0][2] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[0][2]);		
		Fix::virial[5] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[1][2] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[1][2]);		
		
	}

	
	//free memory of derivatives of extended variables
	for(iex=0;iex<nex;iex++){
		for(j=0;j<nmax;j++){
			delete [] exvar_tmp[iex].dQ[j];
		}
		delete [] exvar_tmp[iex].dQ;
	}
}







void FixDafedgen::get_pathcv_values_test()
{
	int i,j;
	int iex,nex;
	int L,mi;
	//double box[3];				//size of simulation box
	double box[6];				//size of simulation box
	double **x = atom->x;			//atom vector, contains only local atoms, owned and ghost
	tagint *tag = atom->tag;
	tagint *moltag = atom->molecule;
	int *type = atom->type;
	DAFED::EXvariables exvar_tmp[2];
	int natoms=int(atom->natoms);
	nex = 2;

	int jnum;
	int *jlist;
	for(iex=0;iex<nex;iex++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				exvar_tmp[iex].dvirial[i][j] = 0.0;
			}
		}
	}
	
	//allocate memory for extended variable derivatives
	for(iex=0;iex<nex;iex++){
		//allocate derivatives
		exvar_tmp[iex].dQ = new double* [nmax];
		for(int j=0;j<nmax;j++){
			exvar_tmp[iex].dQ[j] = new double [3];
		}
	}
	//NNout is used for the trajectory data output. 
	for (i = 0;i<nmax;i++)        //loop over a maximum number of atoms
	{
		for (int outi = 0; outi < parameter.nnout; outi++)    //loop over number of NNoutput
		{                                       //set up values to be zero for NNout 
			molecules[i].NNout[outi] = 0;
			NNout[i][outi] = 0;
		}//end of NNout loop
	}//end of nmax loop	
	// only orthorhombic box works here!
	parameter.triclinic = domain->triclinic;
	box[0] = domain->h[0];
	box[1] = domain->h[1];
	box[2] = domain->h[2];
	if (parameter.triclinic)
	{
		box[3] = domain->h[3];//yz
        	box[4] = domain->h[4];//xz
       		box[5] = domain->h[5];//xy
	}
	else
	{
		box[3] = 0;
        	box[4] = 0;
       		box[5] = 0;	
	}
	//initialize data structure that hold info on lammps neighbour list
	lmpneigh.inum = list->inum;			//same as local num except using  pair_style hybrid. In this case, number of neigh
	lmpneigh.ilist = list->ilist;			//get local id for this specific ato
	lmpneigh.numneigh = list->numneigh;		//number of neighbors
	lmpneigh.firstneigh = list->firstneigh;		//list of neighbors for one specific atom
	lmpneigh.nall = atom->nlocal + atom->nghost;	//number of atoms in the processor(includes ghost)
	lmpneigh.me = me;				//proc number
	lmpneigh.ntot = atom->natoms;			//number of total atoms

	//copy current positions into data structure
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			molecules[i].pos[j] = x[i][j];	//get position for each atom
		}
		molecules[i].moltype = moltag[i]-1;	//get mol type for each atom
		molecules[i].eletype = type[i];		//get element type for each atom
	}
	
//	auto timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
//	cout << "Before all calc " << ctime(&timenow) << endl;

	//get all the neighbor info to molecule class.
	DAFED::get_allNeighbourDistances_lmp_new(molecules,molec,lmpneigh,lmpneigh.nall,box,parameter);
	//comm->forward_comm_fix(this);//Need to communicate here to calculate COvec and NNvec
	comm->forward_comm(this);//Need to communicate here to calculate COvec and NNvec
//	MPI_Barrier(world);
	
//	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
//	cout << "After neighbour " << ctime(&timenow) << endl;
	
	//get symmetry fuction and Neural Network result from here. 
	DAFED::get_vec_symmetry_NN_3_lmpneighbour_new(molecules, molec,lmpneigh,box,parameter,Rskappa,eta,NNvar);
	
//	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
//	cout << "After get_vec_symmetry" << ctime(&timenow) << endl;
	
	//comm->forward_comm_fix(this);//Need to communicate here to calculate NNgrad
	comm->forward_comm(this);//Need to communicate here to calculate NNgrad
//	MPI_Barrier(world);
	
	
	//fwd all derivative info from here.
	parameter.nex = 2;
	DAFED::get_Q_vec_lmpneigh(molecules,molec, lmpneigh,exvar_tmp,parameter,Rskappa,eta,NNvar,Qlocal);
	
//	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
//	cout << "After get Q " << ctime(&timenow) << endl;
//	exit(1);	
	parameter.nex = 1;
	MPI_Barrier(world);

	//sum contributions from all procs. get Q.
	for(i=0;i<nex;i++){
		MPI_Allreduce(&exvar_tmp[i].Qlocal,&exvar_tmp[i].Q,1,MPI_DOUBLE,MPI_SUM,world);
	}
	MPI_Allreduce(&Qlocal,&Qall,parameter.nnout,MPI_DOUBLE,MPI_SUM,world);
	
	
	double QU1[10],QU4[10];		//CHECK THAT THIS IS CONSISTENT WITH NO. OF POINTS
	double DISTQ[10];
	double expD[10];
	double lambda;			//setted up as 200
	double sumN,sumD;		//sum of numerator part in S(Y) and sum of denominator part in S(Y).
	double pathcv;			//S(Y)
	double dpcv_dU1, dpcv_dU4;
	int pathpts;			//No. of points along pathcv
//kokodemo dame
	//set how many points along the pathcv
	pathpts = 10;

	//define images along the path, state 0 is 0.40-U1 and 0.60-U4, last state is 0.80-U1, 0.20-U4 with n images
	for(i=0;i<pathpts;i++){
		//QU1[i] = 0.40 + i*0.04;
		//QU4[i] = 0.55 - i*0.04;
		//looks like U1 is less stable, so try to increase U4 instead.
		//QU1[i] = 0.96 - i*0.07;
		//QU4[i] = 0.14 + i*0.07;
		//QU1[i] = 0.8 - i*0.06;
		//QU4[i] = 0.2 + i*0.06;
		//QU1[i] = 0.40 - i*0.03;
		//QU4[i] = 0.55 + i*0.03;
		QU1[i] = 0.9 - i*0.07;
		QU4[i] = 0.1 + i*0.07;
	}

	lambda = 200.0;

	//exvar_tmp[0] = BCC, exvar_tmp[1] = A15
	//D1 = sqrt((QU1_1-exvar_tmp[0].Q)*(QU1_1-exvar_tmp[0].Q) + (QU4_1-exvar_tmp[1].Q)*(QU4_1-exvar_tmp[1].Q));
	//D2 = sqrt((QU1_2-exvar_tmp[0].Q)*(QU1_2-exvar_tmp[0].Q) + (QU4_2-exvar_tmp[1].Q)*(QU4_2-exvar_tmp[1].Q));
	for(i=0;i<pathpts;i++){
		//DISTQ[i] = sqrt((QU1[i]-exvar_tmp[0].Q)*(QU1[i]-exvar_tmp[0].Q) + (QU4[i]-exvar_tmp[1].Q)*(QU4[i]-exvar_tmp[1].Q));
		DISTQ[i] = (QU1[i]-exvar_tmp[0].Q)*(QU1[i]-exvar_tmp[0].Q) + (QU4[i]-exvar_tmp[1].Q)*(QU4[i]-exvar_tmp[1].Q);	//it is more likely diffQ^2. diffQ=sqrt(diffQU1**2 + diffQU4**2)
		expD[i] = exp(-lambda*DISTQ[i]);		//denominator part(except sum)
	}

	sumN = 0.0;
	sumD = 0.0;
	for(i=0;i<pathpts;i++){
		sumN += i*expD[i];
		sumD += expD[i];
	}
	pathcv = 1.0/(pathpts-1.0) * sumN/sumD;	//S(Y)

	//assign value to cv used in simulation
//	exvar[1].Q = pathcv;		//This is not in tmp!!!
	exvar[0].Q = pathcv;		//This is not in tmp!!!
	//cout << "pathcv = " << pathcv << endl; 

	//now get the derivative as dS/dx = dS/dU1*dU1/dx + dS/dU4*dU4/dx
	//dpcv_dU1 = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QU1_1-exvar_tmp[0].Q)/D1 + expD2*lambda*(QU1_2-exvar_tmp[0].Q)/D2))
	//		/((expD1 + expD2)*(expD1 + expD2)))
	//	+ (expD1*lambda*(QU1_1-exvar_tmp[0].Q)/D1 + 2.0*expD2*lambda*(QU1_2-exvar_tmp[0].Q)/D2)/(expD1 + expD2);

	//dpcv_dU4 = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QU4_1-exvar_tmp[1].Q)/D1 + expD2*lambda*(QU4_2-exvar_tmp[1].Q)/D2))
	//		/((expD1 + expD2)*(expD1 + expD2)))
	//	+ (expD1*lambda*(QU4_1-exvar_tmp[1].Q)/D1 + 2.0*expD2*lambda*(QU4_2-exvar_tmp[1].Q)/D2)/(expD1 + expD2);
	double dsumN_dU1,dsumD_dU1;
	double dsumN_dU4,dsumD_dU4;
	dsumN_dU1 = 0.0;
	dsumD_dU1 = 0.0;
	dsumN_dU4 = 0.0;
	dsumD_dU4 = 0.0;
	for(i=0;i<pathpts;i++){
		dsumN_dU1 += i*lambda*2.0*(QU1[i]-exvar_tmp[0].Q)*expD[i];
		dsumD_dU1 += lambda*2.0*(QU1[i]-exvar_tmp[0].Q)*expD[i];

		dsumN_dU4 += i*lambda*2.0*(QU4[i]-exvar_tmp[1].Q)*expD[i];
		dsumD_dU4 += lambda*2.0*(QU4[i]-exvar_tmp[1].Q)*expD[i];
	}

	dpcv_dU1 = 1.0/(pathpts-1.0)*(dsumN_dU1/sumD - sumN*dsumD_dU1/(sumD*sumD));
	dpcv_dU4 = 1.0/(pathpts-1.0)*(dsumN_dU4/sumD - sumN*dsumD_dU4/(sumD*sumD));

	//for(i=0;i<atom->nlocal;i++){
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			exvar[0].dQ[i][j] = dpcv_dU1 * exvar_tmp[0].dQ[i][j] + dpcv_dU4 * exvar_tmp[1].dQ[i][j];
		}
	}
	//DK:added for NNout.  
	for (int mi = 0; mi < lmpneigh.lmpmol; mi++)  //loop over lmpmol which is a number of the molecule in the proc
	{
		for (int atmi = 0; atmi <parameter.natm; atmi++)     //loop over a number of the atoms in each molecule 
        	{
                	for (int outi = 0; outi < parameter.nnout; outi++){   //loop over number of NNout
                        	molecules[molec[mi].atom_ID[atmi]].NNout[outi] = molecules[molec[mi].atom_ID[0]].NNout[outi];
//				if (mi == 0)
//				{
//					cout << "molecules[molec[mi].atom_ID[atmi]].NNout[outi] = " << molecules[molec[mi].atom_ID[atmi]].NNout[outi] << endl;
//				}
	                }//end of NNout
		}//end of atmi loop
	}//end of lmpmol loop
	//comm->reverse_comm_fix(this);	//get all the ghost info to real atoms. 
	comm->reverse_comm(this);	//get all the ghost info to real atoms. 
	//DK:Trnasfering all the NNout info to the atom_array. 
	for (int ti = 0; ti < atom->nlocal; ti++)   //loop over all owned atoms
	{
		for (i = 0; i < parameter.nnout; i++)         //loop over number of NNout
		{
			NNout[ti][i] =molecules[ti].NNout[i];
		}//end of NNout loop
	}//end of owned atoms loop




	//and also the virial (this is the local sum for each proc, also in lammps these are determined for each proc)
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			exvar[0].dvirial[i][j] = (dpcv_dU1 * exvar_tmp[0].dvirial[i][j]) + (dpcv_dU4 * exvar_tmp[1].dvirial[i][j]);
			//cout << exvar[0].dvirial[i][j] << endl;
		}
	}

	//check if additional restraint is applied (here this is the transition tube)

//**********************************Most likely there is a bug in this part. Be careful****************************	
	if(parameter.use_restraint == 1){
		//cout << "parameter.use_restraint == 1 " << endl;
		lambda = 1000.0;		//larger lambda value
		double dist_pathcv;	//"distance" from pathcv: -1/lambda*Log(sumD)*30
		double diffrest, diffrest_prod;
		double ddist_pathcv_dU1, ddist_pathcv_dU4;

		sumD = 0.0;
		for(i=0;i<pathpts;i++){
			expD[i] = exp(-lambda*DISTQ[i]);
			sumD += expD[i];
		}

		dist_pathcv = -30.0/lambda*log(sumD);


		diffrest = dist_pathcv - parameter.restraint_value;
		diffrest_prod = 1.0;
		for(int iwall=0; iwall < (parameter.restraint_n-1); iwall++){
			diffrest_prod *= diffrest;
		}
		//exvar_tmp[0].frestrain = parameter.restraint_pref * diffrest_prod;	//-dVres/dq
		//exvar_tmp[1].frestrain = parameter.restraint_pref * diffrest_prod;
		exvar[0].frestrain = parameter.restraint_pref * diffrest_prod;		//-dVres(s)/ds  where s=dist_pathcv
		exvar[0].vrestrain = -exvar[0].frestrain * diffrest / parameter.restraint_n;	// Vres(s) 

		//get derivative of dist_pathcv wrt U1 and U4
		dsumD_dU1 = 0.0;
		dsumD_dU4 = 0.0;
		for(i=0;i<pathpts;i++){
			dsumD_dU1 += lambda*2.0*(QU1[i]-exvar_tmp[0].Q)*expD[i];
			dsumD_dU4 += lambda*2.0*(QU4[i]-exvar_tmp[1].Q)*expD[i];
		}

		ddist_pathcv_dU1 = -30.0/(lambda*sumD)*dsumD_dU1;
		ddist_pathcv_dU4 = -30.0/(lambda*sumD)*dsumD_dU4;

		//start: this is just for testing the derivative
		//for(i=0;i<atom->nlocal;i++){
		//	for(j=0;j<3;j++){
		//		exvar[0].dQ[i][j] = -exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dA15*exvar_tmp[1].dQ[i][j]);
		//	}
		//}
		
		//exvar[0].Q = -parameter.restraint_pref/parameter.restraint_n*diffrest_prod*diffrest;

		//end: only for testing derivatives
		
		//add forces to atom
//		cout << "Before adding force to the atom " << endl;
//		DK: CHECK HERE AGAIN. Add forces to atoms crash the code. 
		for(i=0;i<atom->nlocal;i++){
			for(j=0;j<3;j++){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//DK:Most Likely THIS IS THE PROBLEM When use rstraint is 1
				atom->f[i][j] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dU4*exvar_tmp[1].dQ[i][j]);
				//to print out forces due to restrain
				//0-7 are the other forces!
				array[i][8+j] = exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dU4*exvar_tmp[1].dQ[i][j]);
			}
		}
//		cout << "After adding force " << endl;
//DK:commented out for now		//update virial
		Fix::virial[0] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[0][0] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[0][0]);		
		Fix::virial[1] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[1][1] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[1][1]);		
		Fix::virial[2] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[2][2] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[2][2]);		
		Fix::virial[3] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[0][1] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[0][1]);		
		Fix::virial[4] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[0][2] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[0][2]);		
		Fix::virial[5] += exvar[0].frestrain * (ddist_pathcv_dU1*exvar_tmp[0].dvirial[1][2] + ddist_pathcv_dU4*exvar_tmp[1].dvirial[1][2]);		
		
	}

	
	//free memory of derivatives of extended variables
	for(iex=0;iex<nex;iex++){
		for(j=0;j<nmax;j++){
			delete [] exvar_tmp[iex].dQ[j];
		}
		delete [] exvar_tmp[iex].dQ;
	}
}

//------------------------------------------------------------------------
//  calculate the order parameter in two steps, since we need information
//  of neighbours to calculate the derivatives of the Steinhardt parameters
//  get path collective variable as a function of QBCC and QA15
//------------------------------------------------------------------------
/*void FixDafed::get_pathcv_values_Xpoints()
{
	int i,j;
	int iex,nex;
	int L,mi;
	double box[3];				//size of simulation box
	double **x = atom->x;			//atom vector, contains only local atoms, owned and ghost
	tagint *tag = atom->tag;
	DAFED::EXvariables exvar_tmp[2];

	//initialize some stuff for temporary exvar
	nex = 2;
	for(iex=0;iex<nex;iex++){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				exvar_tmp[iex].dvirial[i][j] = 0.0;
			}
		}
	}
	//allocate memory for extended variable derivatives
	for(iex=0;iex<nex;iex++){
		//allocate derivatives
		exvar_tmp[iex].dQ = new double* [nmax];
		for(int j=0;j<nmax;j++){
			exvar_tmp[iex].dQ[j] = new double [3];
		}
	}

	
	// only orthorhombic box works here!
	box[0] = domain->h[0];
	box[1] = domain->h[1];
	box[2] = domain->h[2];

	//initialize data structure that hold info on lammps neighbour list
	lmpneigh.inum = list->inum;
	lmpneigh.ilist = list->ilist;
	lmpneigh.numneigh = list->numneigh;
	lmpneigh.firstneigh = list->firstneigh;
	lmpneigh.nall = atom->nlocal + atom->nghost;
	lmpneigh.me = me;
	lmpneigh.ntot = atom->natoms;

	//copy current positions into data structure
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			molecules[i].pos[j] = x[i][j];
		}
	}
	
	//make small neighbour lists for calculating collective variable 
	DAFED::get_symmetry_NN_3_lmpneighbour(molecules,lmpneigh,box,parameter,sf2Rs,sf2eta,sf3kappa,NNvar);


	//from fcc/orient: copy stuff into nbr structure and communicate
	for(i=0;i<atom->nlocal;i++){
		for(j=0;j<molecules[i].n_neighbors1;j++){
			nbr[i].id[j] = molecules[i].neighbors1[j];
		}
	}

	// communicate to acquire nbr and molecules data for ghost atoms
	comm->forward_comm_fix(this);

	//call second part after communication
	//this is specifically for nex=2 with 0=BCC and 1=A15
	//and                  for nex=1 with 0=A15
	//set parameter.nex = 2 since we want both QBCC and QA15
	parameter.nex = 2;
	DAFED::get_Q_3_lmpneighbour(molecules,lmpneigh,exvar_tmp,parameter,sf2Rs,sf2eta,sf3kappa,NNvar,Qlocal);
	//reset parameter.nex to 1 since only a singel path collective variable
	parameter.nex = 1;
	MPI_Barrier(world);

	//sum contributions from all procs
	for(i=0;i<nex;i++){
		MPI_Allreduce(&exvar_tmp[i].Qlocal,&exvar_tmp[i].Q,1,MPI_DOUBLE,MPI_SUM,world);
	}
	MPI_Allreduce(&Qlocal,&Qall,5,MPI_DOUBLE,MPI_SUM,world);

	
	//here we have now the values of Q-BCC, Q-A15 as well as all dQ-BCC/dx and dQ-A15/dx values
	//double QBCC_1, QBCC_2;		//values in state 1 and 2
	//double QA15_1, QA15_2;
	double QBCC[10],QA15[10];		//CHECK THAT THIS IS CONSISTENT WITH NO. OF POINTS
	//double D1, D2;		//distance D1 and D2 in terms of QBCC and QA15 values;
	double DISTQ[10];
	//double expD1, expD2;
	double expD[10];
	double lambda;
	double sumN,sumD;
	double pathcv;
	double dpcv_dBCC, dpcv_dA15;
	int pathpts;			//No. of points along pathcv

	//set how many points along the pathcv
	pathpts = 10;

	//QBCC_1 = 0.2;		//state 1: 0.2-BCC and 0.6-A15
	//QBCC_2 = 0.6;		//state 2: 0.6-BCC and 0.2-A15
	//QA15_1 = 0.6;
	//QA15_2 = 0.2;
	
	//define images along the path, state 0 is 0.2-BCC and 0.6-A15, last state is 0.65-BCC, 0.15-A15 with n images
	//we have also used 0.2-BCC and 0.56-A15 to 0.65-BCC and 0.11-A15
	//                  0.2-BCC and 0.5-BCC to 0.65-BCC and 0.05-A15
	for(i=0;i<pathpts;i++){
		QBCC[i] = 0.20 + i*0.05;
		QA15[i] = 0.50 - i*0.05;
	}

	lambda = 200.0;

	//exvar_tmp[0] = BCC, exvar_tmp[1] = A15
	//D1 = sqrt((QBCC_1-exvar_tmp[0].Q)*(QBCC_1-exvar_tmp[0].Q) + (QA15_1-exvar_tmp[1].Q)*(QA15_1-exvar_tmp[1].Q));
	//D2 = sqrt((QBCC_2-exvar_tmp[0].Q)*(QBCC_2-exvar_tmp[0].Q) + (QA15_2-exvar_tmp[1].Q)*(QA15_2-exvar_tmp[1].Q));
	for(i=0;i<pathpts;i++){
		//DISTQ[i] = sqrt((QBCC[i]-exvar_tmp[0].Q)*(QBCC[i]-exvar_tmp[0].Q) + (QA15[i]-exvar_tmp[1].Q)*(QA15[i]-exvar_tmp[1].Q));
		DISTQ[i] = (QBCC[i]-exvar_tmp[0].Q)*(QBCC[i]-exvar_tmp[0].Q) + (QA15[i]-exvar_tmp[1].Q)*(QA15[i]-exvar_tmp[1].Q);
		expD[i] = exp(-lambda*DISTQ[i]);
	}

	//expD1 = exp(-lambda*D1);
	//expD2 = exp(-lambda*D2);

	//pathcv = (expD1 + 2.0*expD2)/(expD1 + expD2); 
	sumN = 0.0;
	sumD = 0.0;
	for(i=0;i<pathpts;i++){
		sumN += i*expD[i];
		sumD += expD[i];
	}
	pathcv = 1.0/(pathpts-1.0) * sumN/sumD;

	//assign value to cv used in simulation
	exvar[0].Q = pathcv;

	//now get the derivative as dS/dx = dS/dBCC*dBCC/dx + dS/dA15*dA15/dx
	//dpcv_dBCC = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QBCC_1-exvar_tmp[0].Q)/D1 + expD2*lambda*(QBCC_2-exvar_tmp[0].Q)/D2))
	//		/((expD1 + expD2)*(expD1 + expD2)))
	//	+ (expD1*lambda*(QBCC_1-exvar_tmp[0].Q)/D1 + 2.0*expD2*lambda*(QBCC_2-exvar_tmp[0].Q)/D2)/(expD1 + expD2);

	//dpcv_dA15 = -(((expD1 + 2.0*expD2) * (expD1*lambda*(QA15_1-exvar_tmp[1].Q)/D1 + expD2*lambda*(QA15_2-exvar_tmp[1].Q)/D2))
	//		/((expD1 + expD2)*(expD1 + expD2)))
	//	+ (expD1*lambda*(QA15_1-exvar_tmp[1].Q)/D1 + 2.0*expD2*lambda*(QA15_2-exvar_tmp[1].Q)/D2)/(expD1 + expD2);
	double dsumN_dBCC,dsumD_dBCC;
	double dsumN_dA15,dsumD_dA15;
	dsumN_dBCC = 0.0;
	dsumD_dBCC = 0.0;
	dsumN_dA15 = 0.0;
	dsumD_dA15 = 0.0;
	for(i=0;i<pathpts;i++){
		dsumN_dBCC += i*lambda*2.0*(QBCC[i]-exvar_tmp[0].Q)*expD[i];
		dsumD_dBCC += lambda*2.0*(QBCC[i]-exvar_tmp[0].Q)*expD[i];

		dsumN_dA15 += i*lambda*2.0*(QA15[i]-exvar_tmp[1].Q)*expD[i];
		dsumD_dA15 += lambda*2.0*(QA15[i]-exvar_tmp[1].Q)*expD[i];
	}

	dpcv_dBCC = 1.0/(pathpts-1.0)*(dsumN_dBCC/sumD - sumN*dsumD_dBCC/(sumD*sumD));
	dpcv_dA15 = 1.0/(pathpts-1.0)*(dsumN_dA15/sumD - sumN*dsumD_dA15/(sumD*sumD));

	for(i=0;i<atom->nlocal;i++){
		for(j=0;j<3;j++){
			exvar[0].dQ[i][j] = dpcv_dBCC * exvar_tmp[0].dQ[i][j] + dpcv_dA15 * exvar_tmp[1].dQ[i][j];
		}
	}

	//and also the virial (this is the local sum for each proc, also in lammps these are determined for each proc)
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			exvar[0].dvirial[i][j] = (dpcv_dBCC * exvar_tmp[0].dvirial[i][j]) + (dpcv_dA15 * exvar_tmp[1].dvirial[i][j]);
		}
	}

	//check if additional restraint is applied (here this is the transition tube)
	if(parameter.use_restraint == 1){
		lambda = 1000.0;		//larger lambda value
		double dist_pathcv;	//"distance" from pathcv: -1/lambda*Log(sumD)*30
		double diffrest, diffrest_prod;
		double ddist_pathcv_dBCC, ddist_pathcv_dA15;

		sumD = 0.0;
		for(i=0;i<pathpts;i++){
			expD[i] = exp(-lambda*DISTQ[i]);
			sumD += expD[i];
		}

		dist_pathcv = -30.0/lambda*log(sumD);


		diffrest = dist_pathcv - parameter.restraint_value;
		diffrest_prod = 1.0;
		for(int iwall=0; iwall < (parameter.restraint_n-1); iwall++){
			diffrest_prod *= diffrest;
		}
		//exvar_tmp[0].frestrain = parameter.restraint_pref * diffrest_prod;	//-dVres/dq
		//exvar_tmp[1].frestrain = parameter.restraint_pref * diffrest_prod;
		exvar[0].frestrain = parameter.restraint_pref * diffrest_prod;		//-dVres(s)/ds  where s=dist_pathcv
		exvar[0].vrestrain = -exvar[0].frestrain * diffrest / parameter.restraint_n;	// Vres(s) 

		//get derivative of dist_pathcv wrt BCC and A15
		dsumD_dBCC = 0.0;
		dsumD_dA15 = 0.0;
		for(i=0;i<pathpts;i++){
			dsumD_dBCC += lambda*2.0*(QBCC[i]-exvar_tmp[0].Q)*expD[i];
			dsumD_dA15 += lambda*2.0*(QA15[i]-exvar_tmp[1].Q)*expD[i];
		}

		ddist_pathcv_dBCC = -30.0/(lambda*sumD)*dsumD_dBCC;
		ddist_pathcv_dA15 = -30.0/(lambda*sumD)*dsumD_dA15;

		//start: this is just for testing the derivative
		//for(i=0;i<atom->nlocal;i++){
		//	for(j=0;j<3;j++){
		//		exvar[0].dQ[i][j] = -exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dA15*exvar_tmp[1].dQ[i][j]);
		//	}
		//}
		
		//exvar[0].Q = -parameter.restraint_pref/parameter.restraint_n*diffrest_prod*diffrest;

		//end: only for testing derivatives
		
		//add forces to atom
		for(i=0;i<atom->nlocal;i++){
			for(j=0;j<3;j++){
				atom->f[i][j] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dA15*exvar_tmp[1].dQ[i][j]);
				//to print out forces due to restrain
				//0-7 are the other forces!
				array[i][8+j] = exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dQ[i][j] + ddist_pathcv_dA15*exvar_tmp[1].dQ[i][j]);
			}
		}
		//update virial
		Fix::virial[0] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dvirial[0][0] + ddist_pathcv_dA15*exvar_tmp[1].dvirial[0][0]);		
		Fix::virial[1] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dvirial[1][1] + ddist_pathcv_dA15*exvar_tmp[1].dvirial[1][1]);		
		Fix::virial[2] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dvirial[2][2] + ddist_pathcv_dA15*exvar_tmp[1].dvirial[2][2]);		
		Fix::virial[3] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dvirial[0][1] + ddist_pathcv_dA15*exvar_tmp[1].dvirial[0][1]);		
		Fix::virial[4] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dvirial[0][2] + ddist_pathcv_dA15*exvar_tmp[1].dvirial[0][2]);		
		Fix::virial[5] += exvar[0].frestrain * (ddist_pathcv_dBCC*exvar_tmp[0].dvirial[1][2] + ddist_pathcv_dA15*exvar_tmp[1].dvirial[1][2]);		
		
	}

	
	//free memory of derivatives of extended variables
	for(iex=0;iex<nex;iex++){
		for(j=0;j<nmax;j++){
			delete [] exvar_tmp[iex].dQ[j];
		}
		delete [] exvar_tmp[iex].dQ;
	}


}
*/



//------------------------------------------------------------------------
//  calculate the order parameter in two steps, since we need information
//  of neighbours to calculate the derivatives of the Steinhardt parameters
//------------------------------------------------------------------------
/*void FixDafed::get_Q6_values()
{
	int i,j;
	int mi;
	double box[3];				//size of simulation box
	double **x = atom->x;			//atom vector, contains only local atoms, owned and ghost
	tagint *tag = atom->tag;

	double realmi[13],imgmi[13];
	double realmi_local[13],imgmi_local[13];
	//double sumQ6;
	int nlocal = int(atom->nlocal);
	int ntotal = int(atom->natoms);

	
	// only orthorhombic box works here!
	box[0] = domain->h[0];
	box[1] = domain->h[1];
	box[2] = domain->h[2];

	//initialize data structure that hold info on lammps neighbour list
	lmpneigh.inum = list->inum;
	lmpneigh.ilist = list->ilist;
	lmpneigh.numneigh = list->numneigh;
	lmpneigh.firstneigh = list->firstneigh;
	lmpneigh.nall = atom->nlocal + atom->nghost;
	lmpneigh.me = me;
	lmpneigh.ntot = atom->natoms;

	//cout<<me<<":  inum = "<<lmpneigh.inum<<"  nlocal = "<<atom->nlocal<<endl;

	//copy current positions into data structure
	for(i=0;i<lmpneigh.nall;i++){
		for(j=0;j<3;j++){
			molecules[i].pos[j] = x[i][j];
		}
	}

	//determine distances and cutoff functions
	DAFED::get_allNeighbourDistances_cutoff_lmpneighbour(molecules,lmpneigh,box,parameter);

	//get realmi and imgmi values for atoms on this proc
	//DAFED::calculate_globalQ6_lmpneighbour_part01(molecules,nlocal,ntotal,realmi_local,imgmi_local);
	DAFED::calculate_globalQ6_virial_lmpneighbour_part01(molecules,nlocal,ntotal,realmi_local,imgmi_local);
	MPI_Allreduce(&realmi_local,&realmi,13,MPI_DOUBLE,MPI_SUM,world);
	MPI_Allreduce(&imgmi_local,&imgmi,13,MPI_DOUBLE,MPI_SUM,world);

	
	//DAFED::calculate_globalQ6_lmpneighbour_part02(molecules,nlocal,ntotal,realmi,imgmi,exvar);
	DAFED::calculate_globalQ6_virial_lmpneighbour_part02(molecules,nlocal,ntotal,realmi,imgmi,exvar);
/
}*/

// ------------------------------------------------------------------------
// Get numerical derivatives and compare to analytical ones
// ------------------------------------------------------------------------
// This is for path CV
void FixDafedgen::get_numerical_derivatives()
{
	
	//get numerical derivatives
	double delta=1e-4;		//displacement for numerical derivative		1e-4
	double numer;			//value of numerical derivative
	double **xorig;			//save original atom positions
	double **dQorig;		//save analytical derivatives
	double Q2P,QP,QM,Q2M;		//Q values with displaced atoms
	int myiex;			//select extended variable
	int nlocal = int(atom->nlocal);	//no of owned atoms on this proc
	int natoms = int(atom->natoms);	//total no of atoms in simulation	
	tagint *tag = atom->tag;	//global atom ids
	double **x = atom->x;		//atom position vector
	int i,j;
	int *type = atom->type;

	myiex = 0;			//0: will calc pathcv. For pathcv, we had only one component. For U1
	//myiex = 1;			//0: will calc pathcv. For pathcv, we had only one component. For liq
	xorig = new double* [natoms];
	dQorig = new double* [natoms];

	for(i=0;i<nlocal;i++){		//lmpneigh.nall original --> nlocal
		xorig[i] = new double[3];
		dQorig[i] = new double[3];
	}
	for(i=0;i<nlocal; i++){			//lets try with natoms, original was nlocal
		for(j=0;j<3;j++){
			xorig[i][j] = atom->x[i][j];
			dQorig[i][j] = exvar[myiex].dQ[i][j];
		}
	}
//////////////////############################
//	myiex = 1;			//0: will calc pathcv. For pathcv, we had only one component. For U1
//////////////////############################
	ofstream f_grad;		//output file for numerical/analytical derivatives
	ostringstream ss;
	string fgradname;
	ss<<me;
	fgradname = "gradient_qlm-"+ss.str()+".txt";	//open a separate file for each proc
	f_grad.open(fgradname.c_str(),ofstream::out);
	f_grad<<"#Gradients from proc "<<me<<endl;

	int nall = atom->nlocal+atom->nghost;		//no of owend + ghost atoms on this proc
		
	for(int ti=1;ti<natoms;ti++){		//loop over all atom ids			//should be for it =1;ti<natomsnatoms, but for testing, ti <10
//THIS IS FOR TESTING
//	for(int ti=720;ti<natoms;ti++){		//loop over SOME atom ids			//should be for it =1;ti<natomsnatoms, but for testing, ti <10
//	for(int ti=360;ti<natoms;ti++){		//loop over SOME atom ids			//should be for it =1;ti<natomsnatoms, but for testing, ti <10
//		if (!(ti >2452 && ti < 2460)) continue;	//This was for testing
//		if (!(ti >4100 && ti < 4300)) continue;	//This was for testing
//		if(type[tag[ti]]==1)continue;	//IF you wont to check for hydrogen, remove this line. If not, this will save almost half of the simulation time for urea.
		//for(int ti=1;ti<100;ti++){		//loop over all atom ids			//should be for it =1;ti<natomsnatoms, but for testing, ti <10
		for(j=0;j<3;j++){		//loop over x,y,z component
			//get Q2P
			for(int iproc=0;iproc<nprocs;iproc++){		//find which proc has this atom and modify component
				if(iproc == me){
					for(i=0;i<nlocal;i++){
						if(tag[i] == ti){
							x[i][j] = xorig[i][j] + 2.0*delta;
						}
					}
				}
				MPI_Barrier(world);
			}
			comm->forward_comm();
			MPI_Barrier(world);
			//now call function to calculate Q
			if(parameter.choose_exvar == 0) cout << "ERROR";/*get_exvar_values();*/
			else if(parameter.choose_exvar == 1) get_pathcv_values_test();
			Q2P = exvar[myiex].Q;		//this should be fine on all procs

			//get QP
			for(int iproc=0;iproc<nprocs;iproc++){		//repeat for next displacement
				if(iproc == me){
					for(i=0;i<nlocal;i++){
						if(tag[i] == ti){
							x[i][j] = xorig[i][j] + 1.0*delta;
						}
					}
				}
				MPI_Barrier(world);
			}
			comm->forward_comm();
			MPI_Barrier(world);
			if(parameter.choose_exvar == 0) cout << "ERROR";//get_exvar_values();
			else if(parameter.choose_exvar == 1) get_pathcv_values_test();
			QP = exvar[myiex].Q;
			//get QM
			for(int iproc=0;iproc<nprocs;iproc++){
				if(iproc == me){
					for(i=0;i<nlocal;i++){
						if(tag[i] == ti){
							x[i][j] = xorig[i][j] - 1.0*delta;
						}
					}
				}
				MPI_Barrier(world);
			}
			comm->forward_comm();
			MPI_Barrier(world);
			if(parameter.choose_exvar == 0) cout << "ERROR";//get_exvar_values();
			else if(parameter.choose_exvar == 1) get_pathcv_values_test();
			QM = exvar[myiex].Q;
			//get Q2M	
			for(int iproc=0;iproc<nprocs;iproc++){
				if(iproc == me){
					for(i=0;i<nlocal;i++){
						if(tag[i] == ti){
							x[i][j] = xorig[i][j] - 2.0*delta;
						}
					}
				}
				MPI_Barrier(world);
			}
			comm->forward_comm();
			MPI_Barrier(world);
			if(parameter.choose_exvar == 0) cout << "ERROR";//get_exvar_values();get_exvar_values();
			else if(parameter.choose_exvar == 1) get_pathcv_values_test();
			Q2M = exvar[myiex].Q;	

			//reset to original position
			for(int iproc=0;iproc<nprocs;iproc++){
				if(iproc == me){
					for(i=0;i<nlocal;i++){
						if(tag[i] == ti){
							x[i][j] = xorig[i][j];		//reset position 
						}
					}
				}
				MPI_Barrier(world);
			}
			comm->forward_comm();

			numer = 1.0/(12.0*delta)*(-Q2P+8.0*QP-8.0*QM+Q2M);	//calculate numerical derivative
			for(i=0;i<nlocal;i++){					//print analytical/numerical derivative and compare
//				if(type[i]==1)continue;
				if(tag[i] == ti){
					f_grad<<"atom "<<i<<" id is " << tag[i]<< " type: " << type[i] << "  comp "<<j<<"  ";
					//f_grad<<"pos  "<<x[i][j]<<"  ";
					f_grad<<dQorig[i][j]<<"  "<<numer<<"  ";
					if(fabs(dQorig[i][j]-numer)>1e-10){
						f_grad<<"  1";
					}
					else{
						f_grad<<"  0";
					}
					if(fabs(dQorig[i][j]-numer)>1e-6){
						f_grad<<"  1";
					}
					else{
						f_grad<<"  0";
					}
					if(fabs((dQorig[i][j]-numer)/dQorig[i][j]) > 1e-4){	//relative error 0.01%
						f_grad<<"  1";
					}
					else{
						f_grad<<"  0";
					}
					f_grad<<endl;
				}
			}
		} // end loop over j (x,y,z)
		//exit(0);
	} // end loop over atom ids
	f_grad.close();			//close output file

	for(i=0;i<nlocal;i++){		//delete allocated memory
		delete [] xorig[i];
		delete [] dQorig[i];
	}
	delete [] xorig;
	delete [] dQorig;
	//end numerical derivatives
}



// ------------------------------------------------------------------------
//  from fcc/orient: memory usage of local atom-based arrays
//------------------------------------------------------------------------- 
double FixDafedgen::memory_usage()
{
	// !!THIS IS NOT CORRECT, NEED TO ADD SIZE OF MOLECULES...!!
	double bytes = nmax * sizeof(Nbr);
	//bytes += 2*nmax * sizeof(double);  //that's for 'order' double array
	return bytes;
}

// ---------------------------------------------------------------------- 
//  from fcc/orient: pack data for forward communication
// ---------------------------------------------------------------------- 
//CHECK IF I AM MISSING SOMETHING HERE <-Sample had Steindhart, but I dont use that
int FixDafedgen::pack_forward_comm(int n, int *list, double *buf,int pbc_flag, int *pbc)
{
	
  	int i,j,k,num;
	int L, mi;
  	tagint id;
  	tagint *tag = atom->tag;
  	int nlocal = atom->nlocal;
  	int m = 0;

	for (i = 0; i < n; i++) {
		k = list[i];
		for(j=0;j<3;j++){
			buf[m++] = molecules[k].COvec[j];
			buf[m++] = molecules[k].NNvec[j];
		}
		for(j=0;j<NNvar.output*NNvar.input;j++){
			buf[m++] = molecules[k].NNgrad[j];
		}
	}
	return m;
}
// ----------------------------------------------------------------------
//  from fcc/orient:  unpack data from forward communication
// ----------------------------------------------------------------------
void FixDafedgen::unpack_forward_comm(int n, int first, double *buf)
{
  	int i,j,k,num;
	int L,mi;
  	int last = first + n;
  	int m = 0;

	for (i = first; i < last; i++) {
		for(j=0;j<3;j++){
			molecules[i].COvec[j] = buf[m++];
			molecules[i].NNvec[j] = buf[m++];
		}
		for(j=0;j<NNvar.output*NNvar.input;j++){
			molecules[i].NNgrad[j] = buf[m++];
                }
	}
}

int FixDafedgen::pack_reverse_comm(int n, int first, double *buf)
{
	int i,k,num,m,last;

	m = 0;
	last = first + n;
	//for (i = first; i < last; i++) {
	for (i = first; i < last; i++) {
		for (int k = 0; k< 3;k++){
			buf[m++] = exvar[0].dQ[i][k];
			//buf[m++] = exvar[1].dQ[i][k];
		}
		for (int k = 0; k< parameter.nnout;k++){
			buf[m++] = molecules[i].NNout[k];
		}
	}
	return m;
}

void FixDafedgen::unpack_reverse_comm(int n, int *list, double *buf)
{
	int i,j,m;

	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];	
		for (int k = 0; k< 3;k++){
			exvar[0].dQ[j][k] += buf[m++];
			//exvar[1].dQ[j][k] += buf[m++];
		}
		for (int k = 0; k< parameter.nnout;k++){
			molecules[j].NNout[k] += buf[m++];
		}
	}
}

// ---------------------------------------------------------------------- 
// ---------------------------------------------------------------------- 
// ------------------------------------------------------------------------
// Setup histograms etc for additional metadynamics bias potential / UFED
// ------------------------------------------------------------------------
void FixDafedgen::setup_ufedbias()
{
	int nex = parameter.nex;
	int *index;			//get index for each exvar
	index = new int[nex];
	
	if(parameter.useMetadyn == 1){
		of<<"\nPerforming UFED simulation with additional bias:\n";
	}
	else if (parameter.useMetadyn == 2){
		of<<"\nPerforming Metadynamics simulation with additional bias:\n";
	}
	//else error->all(FLERR,"Parameter use_metadyn must be either 0 = no, 1 = UFED, or 2 = Metadyanmics bias!");
	//else error->universe_one(FLERR,"Parameter use_metadyn must be either 0 = no, 1 = UFED, or 2 = Metadyanmics bias!");
	else{
		cout<<"\n\n!!ERROR!! Parameter use_metadyn must be either 0 = no, 1 = UFED, or 2 = Metadyanmics bias!\n\n";
		exit(1);
	}

	Nlattot = 1;
	Nvertex = 1;
	//setup histogram parameters
	for(int iex = 0; iex<nex; iex++){
		exvar[iex].histo.nlat = exvar[iex].histo.nbin + 1;
		Nlattot *= exvar[iex].histo.nlat;
		exvar[iex].histo.binwidth = (exvar[iex].histo.max - exvar[iex].histo.min)/((double)exvar[iex].histo.nbin);
		if(exvar[iex].histo.binwidth <= 0){
			cout<<"\n\n\t!!ERROR!! Width of bin not well defined for exvar "<<iex<<": "<<exvar[iex].histo.binwidth<<endl;;
			cout<<"histo.min = "<<exvar[iex].histo.min<<"   histo.max = "<<exvar[iex].histo.max<<endl;
			cout<<"Exiting programme...\n";
			exit(1);
		}
		of<<"exvar "<<iex<<": nbin = "<<exvar[iex].histo.nbin
			<<"   ngrid = "<<exvar[iex].histo.nlat
			<<"   binwidth = "<<exvar[iex].histo.binwidth
			<<"   min = "<<exvar[iex].histo.min
			<<"   max = "<<exvar[iex].histo.max
			<<endl;	
		//dimension of neighbour lattice
		Nvertex *= 2;
	}
	of<<"Total no. of grid point = "<<Nlattot<<endl;
	//get factors for converting to/from 1D array
	int chunk;
	chunk = Nlattot;
	for(int iex = 0; iex<nex; iex++){
		chunk /= exvar[iex].histo.nlat;
		exvar[iex].histo.NconvDim = chunk;
	}
	//allocate array for storing bias potential
	bhist = new DAFED::Bhist[Nlattot];
	for(int i=0; i<Nlattot; i++){
		bhist[i].x = new double[nex];
	}	
	//initialize array
	for(int i=0; i<Nlattot; i++){
		bhist[i].pot = 0.0;
		get_ind_rev(index,i);
		for(int iex=0; iex<nex; iex++){
			bhist[i].x[iex] = exvar[iex].histo.min+index[iex]*exvar[iex].histo.binwidth;
		}
	}
	// initialize neighbourhood lattice for interpolation
	// in 2D: f(x,y) = 1/[(x1-x0)*(y1-y0)] * [f(x0,y0)*(x1-x)*(y1-y) + f(x1,y0)*(x-x0)*(y1-y) + f(x0,y1)*(x1-x)*(y-y0) + f(x1,y1)*(x-x0)*(y-y0)]
	//      2:(0,1) -------------- 3:(1,1)
	//              |            |
	//              |            |
  	//              |            |
 	//              |            |
	//              |            |
	//   	0:(0,0) -------------- 1:(1,0)
	//
	// allocate nblat
	nblat = new DAFED::Nblat[Nvertex];
	for(int i=0; i<Nvertex; i++){
		nblat[i].diff = new int[nex];
	}
	int ind, ind_s, sum_exp;
	for(int i=0; i<Nvertex; i++){
		ind = i;
		sum_exp = 0;
		for(int iex=0; iex<nex; iex++){
			ind_s = (ind>>1)<<1;		// bitwise shift to get the right 0/1 indices
			nblat[i].diff[iex] = ind-ind_s;
			ind = ind>>1;			// bitwise shift
			sum_exp += nblat[i].diff[iex];
		}
		if(sum_exp%2==0) nblat[i].sign = 1.0;
		else nblat[i].sign = -1.0;
	}
	//print bias potential properties
	of<<"bias updated every "<<parameter.gauss_freq<<" steps";
	//JR: bias-wall start
	//of<<" AND whenever system is above/below a WALL!";
	//JR: bias-wall end
	of<<endl;
	of<<"Gaussian height = "<<parameter.gauss_h<<"   Gaussian width = "<<parameter.gauss_sigma<<endl;
	of<<"Bias potential printed every "<<parameter.bias_stride<<" steps"<<endl;
	of<<endl;
	
	// set initial bias force = 0 (later also add to read in a particular bias)
	for(int iex = 0; iex<nex; iex++){
		exvar[iex].fbias = 0.0;
	}

	//print no of grid points and min/mas in header of biaspot file
	of_biaspot<<"#No. of grid points = "<<Nlattot;
	for(int iex = 0; iex<nex; iex++){
		of_biaspot<<setprecision(6)<<scientific<<"  min"<<iex<<" = "<<bhist[0].x[iex]<<"  max"<<iex<<" = "<<bhist[Nlattot-1].x[iex];
	}
	of_biaspot<<endl;

	delete [] index;

}

// ------------------------------------------------------------------------
// Update additional metadynamics bias potential / UFED
// ------------------------------------------------------------------------
void FixDafedgen::update_ufedbias()
{
	int nex = parameter.nex;
	double s,x;
	double diff, sumdiff;
	double h = parameter.gauss_h;
	double sigma2 = parameter.gauss_sigma*parameter.gauss_sigma;

	//try to change potential to avoid trapping on upper value of s
	//!!!here only a quick solution for 1D with fixed values!!!  -- NOT generally applicable!!
	//double smax = 0.972;	//this value should be variable


	//if(exvar[0].x < smax){			//only for the first cv!
		for(int i=0; i<Nlattot; i++){
			sumdiff = 0.0;
			for(int iex=0; iex<nex; iex++){
				x = bhist[i].x[iex];
				s = exvar[iex].x;
				diff = fabs(s-x);
				sumdiff += diff*diff;
			}
			bhist[i].pot += h*exp(-sumdiff/(2.0*sigma2)); 
		}
	//}
	//else{				//modified bias to avoid trapping ONLY 1D at the moment!!!
	//	s = smax;
	//	for(int i=0; i<Nlattot; i++){
	//		x = bhist[i].x[0];
	//		diff = fabs(s-x);
	//		sumdiff = diff*diff;

	//		bhist[i].pot += (h-h*exp(-sumdiff/(2.0*sigma2)))*(0.5*tanh((x-s)/sigma2)+0.5);		//inverted gaussion*tanh to be one-sided
	//	}
	//}



}

// ------------------------------------------------------------------------
// Read in additional metadynamics bias potential / UFED
// ------------------------------------------------------------------------
void FixDafedgen::read_ufedbias()
{
	ifstream if_readbias;
	int nex = parameter.nex;
	int iex,i;
	int npoints;
	double *xmin, *xmax;
	double epsi = 1e-8;
	char dummy[2];

	xmin = new double [nex];
	xmax = new double [nex];
	
	if_readbias.open(parameter.bias_readfile.c_str(),ifstream::in);
	if(if_readbias.is_open()){
		//first the # sign
		if_readbias >> dummy[0];
		//read in no. of grid points
		if_readbias >> npoints;
		if(npoints != Nlattot){
			cerr<<"\n\n\t!!ERROR!! No. of grid points in bias histogram not consistent between setup and read-in value:"<<endl;
			cerr<<"\tsetup Nlattot = "<<Nlattot<<"  read-in points = "<<npoints<<endl;
			cerr<<"\tExiting programme..."<<endl<<endl;;
			exit(1);
		}

		for(iex=0;iex<nex;iex++){
			if_readbias >> xmin[iex];
			if_readbias >> xmax[iex];
		}

		for(iex=0;iex<nex;iex++){
			if(fabs(xmin[iex]- bhist[0].x[iex]) > fabs(xmin[iex])*epsi || fabs(xmax[iex]-bhist[Nlattot-1].x[iex]) > fabs(xmax[iex])*epsi){
				cerr<<"\n\n\t!!ERROR!! Min/max values of bias histogram not consistent between setup and read-in values:"<<endl;
				//for (int k=0;k<nex;k++){
				        cerr<<setprecision(6)<<scientific;
					cerr<<"\tsetup  : xmin"<<iex<<" = "<<bhist[0].x[iex]<<"  xmax"<<iex<<" = "<<bhist[Nlattot-1].x[iex]<<endl;
					cerr<<"\tread-in: xmin"<<iex<<" = "<<xmin[iex]<<"  xmax"<<iex<<" = "<<xmax[iex]<<endl;
				//}
				cerr<<"\tExiting programme..."<<endl<<endl;;
				exit(1);
			}
		}
		for(i=0;i<Nlattot;i++){
			for(iex=0;iex<nex;iex++){
				if_readbias >> bhist[i].x[iex];
			}
			if_readbias >> bhist[i].pot;
		}


	}
	else{
		cerr<<"\n\n\t!!ERROR!! Cannot open file '"<<parameter.bias_readfile<<"' to read in bias"<<endl;
		cerr<<"\tExiting programme..."<<endl<<endl;;
		exit(1);
	}

	delete [] xmin;
	delete [] xmax;
	
}


// ------------------------------------------------------------------------
// Determine force from recorded bias potential
// ------------------------------------------------------------------------
void FixDafedgen::get_bias_force()
{
	
	int iex,nex;
	int *index;	// histogram index for each exvar
	int *index_tmp, *index_c;
	int ind, ind_c;
	double pref;
	double pot_v,x_c;
	double prod;

	nex = parameter.nex;
	index = new int[nex];
	index_tmp = new int[nex];
	index_c = new int[nex];
	pref = 1.0;

	for(iex=0;iex<nex;iex++){
		index[iex] = (int)((exvar[iex].x - exvar[iex].histo.min)/exvar[iex].histo.binwidth);  //always lower bound of bin
		if(index[iex] < 0 || index[iex] >= (exvar[iex].histo.nlat-1)){
			of<<"\n\n\t!!ERROR!! Histogram boundaries for bias potential not correct\n";
			of<<"histo index of exvar s"<<iex<<" = "<<index[iex]<<endl;
			of<<"histomin = "<<exvar[iex].histo.min<<"   histomax = "<<exvar[iex].histo.max<<endl;
			of<<"binwidth = "<<exvar[iex].histo.binwidth<<"   no. of grid point = "<<exvar[iex].histo.nlat<<endl;
			of<<"exvar value = "<<exvar[iex].x<<endl;
			of<<endl;
			cout<<"\n\n\tHistogram boundaries for bias potential not correct, check log files and input\n";
			cout<<"\tExiting programme...\n\n";
			exit(1);
		}
		pref *= exvar[iex].histo.binwidth; 	//\prod (x1-x0), only for equal bin width
		exvar[iex].fbias = 0.0;
	}
	//loop over no. of vertices for interpolation
	for(int i=0;i<Nvertex;i++){
		for(iex=0;iex<nex;iex++){
			index_tmp[iex] = index[iex] + nblat[i].diff[iex];
			index_c[iex] = index[iex] + nblat[Nvertex-1-i].diff[iex];
			//cout<<"vertex "<<i<<", index_tmp["<<iex<<"] = "<<index_tmp[iex]<<", index_c["<<iex<<"] = "<<index_c[iex]<<endl;
		}
		//convert to 1D index
		ind = get_ind(index_tmp);		//value at vertex
		ind_c = get_ind(index_c);		//value at 'opposite' vertex
		//cout<<"vertex "<<i<<", ind = "<<ind<<", ind_c = "<<ind_c<<endl;
		if(ind<0 || ind>=Nlattot || ind_c<0 || ind_c>=Nlattot){
			of<<"\n\n\t!!ERROR!! Histogram indices for bias potential not correct\n";
			of<<"vertex = "<<i<<", ind = "<<ind<<", ind_c = "<<ind_c<<endl;
			of<<endl;
			cout<<"\n\n\tHistogram indices for bias potential not correct, check log files and input\n";
			cout<<"\tExiting programme...\n\n";
			exit(1);
		}
		//compute force fbias = -dVbias/dx
		pot_v = bhist[ind].pot;		//value of bias potential at vertex i
		for(iex=0;iex<nex;iex++){
			prod = nblat[i].sign;
			for(int jex=0;jex<nex;jex++){
				if(iex != jex){
					prod *= (bhist[ind_c].x[jex] - exvar[jex].x);	// (x1 - x), resp -(x0 - x)
				}
			}
			exvar[iex].fbias += (pot_v*prod/pref);	//sum over interpolation vertices
		}



	}

	delete [] index;
	delete [] index_tmp;
	delete [] index_c;
	


}

// ------------------------------------------------------------------------
// fixed bias force that changes sign
// ------------------------------------------------------------------------
void FixDafedgen::get_bias_fixedforce()
{
	//int iex, nex;
	double smin,smax;
	double const_force;

	//nex = parameter.nex;
	smin = 0.1;
	smax = 0.9;
	const_force = 300.0;		//this is a bit arbitrary, let's see

	if (exvar[0].x < smin){
		biasforce_sign = 1.0;
	}
	if(exvar[0].x > smax){
		biasforce_sign = -1.0;
	}

	exvar[0].fbias = const_force*biasforce_sign;

}
// ------------------------------------------------------------------------
// convert 1D array into index for each exvar
// ------------------------------------------------------------------------
void FixDafedgen::get_ind_rev(int *index, int i){
	int res, i_dim;
	res = i;
	for(i_dim=0; i_dim<parameter.nex; i_dim++){
		index[i_dim] = (int)(res/exvar[i_dim].histo.NconvDim);
		res = res%exvar[i_dim].histo.NconvDim;
	}
}


// ------------------------------------------------------------------------
// convert n-dim index for n exvar into 1D index
// ------------------------------------------------------------------------
int FixDafedgen::get_ind(int *index){
	int i_dim;
	int ind_1D;

	ind_1D = 0;
	for(i_dim=0; i_dim<parameter.nex;i_dim++){
		ind_1D += index[i_dim]*exvar[i_dim].histo.NconvDim;
	}
	return ind_1D;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixDafedgen::grow_arrays(int nmax)
{
  memory->grow(array,nmax,nvalues,"fix_ave/atom:array");/**/
  //array_atom = array;/**/
  //if (array) vector_atom = array[0];
  //else vector_atom = NULL;
}


/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixDafedgen::copy_arrays(int i, int j, int /*delflag*/)
{
  for (int m = 0; m < nvalues; m++)
    array[j][m] = array[i][m];
}


/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixDafedgen::pack_exchange(int i, double *buf)
{
  for (int m = 0; m < nvalues; m++) buf[m] = array[i][m];
  return nvalues;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixDafedgen::unpack_exchange(int nlocal, double *buf)
{
  for (int m = 0; m < nvalues; m++) array[nlocal][m] = buf[m];
  return nvalues;
}

