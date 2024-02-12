/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_G2G3_gen.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
//#include "compute_chunk_atom.h"
#include "compute_vector_chunk.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "utils.h"
#include "pair.h"

//JR: additional header files
#include <iostream>
#include <iomanip>
//for neighborlists
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
//for communication with ghost atoms
#include "comm.h"

//some math functions for vectors
#include "math_extra.h"


using namespace LAMMPS_NS;
//JR: add another namespace
using namespace std;
using namespace MathExtra;

enum{ONCE,NFREQ,EVERY};

// ---------------------------------------------------------------------- 
// constructor
// ---------------------------------------------------------------------- 
ComputeG2G3gen::ComputeG2G3gen(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idvecchunk1(NULL),
  idvecchunk2(NULL),
  //DK added
  of_g2g3func(NULL),
  me(0),
  nprocs(0),
  ichunk(NULL)
{
	if (narg != 8) error->all(FLERR,"Illegal compute histo/vecchunk command");

	//This is necessary to properly communicate stuff:
  	if (!atom->tag_enable) error->all(FLERR,"compute histo/vecchunk requires atom tags");

	//check if this is run on more than one core
	MPI_Comm_rank(world,&me);
  	MPI_Comm_size(world,&nprocs);
	if(nprocs > 1){
		cout<<"Hello, I'm running compute G2G3gen/vecchunk in parallel with nprocs : "<<nprocs<<endl;
		cout<<"This is proc no "<<me
			<<" nlocal = "<<atom->nlocal
			<<endl;
	}


	// fix the size of the histogram here
	// for now nbin = 201 for -1 <= histo <= 1
	RES = 100.0;
	double m_max, m_min;
	m_max = 1.01;
	m_min = -1.01;
	nbin = int(ceil((m_max - m_min)*RES)) + 1;
	LOW = int(round(m_min*RES));

	//in this compute the 'compute_array()' function is called, values storred in variabl 'array'
	array_flag = 1;	
	size_array_cols = 2;
	size_array_rows = nbin;

	//for communication with ghost atoms
	peratom_flag = 1;
	comm_forward = 1; 	//how many data per atom communicated (here just cluster id)

	//allocate memory for histogram
	memory->create(hist,1,nbin,"histovecchunk:hist");
	memory->create(histall,1,nbin,"histovecchunk:histall");
	memory->create(array,nbin,2,"histvecchunk:array");


	// ID of compute vector/chunk
	
	int n1 = strlen(arg[4]) + 1;
	int n2 = strlen(arg[5]) + 1;   //It should be same as n1 right?
	idvecchunk1 = new char[n1];
	idvecchunk2 = new char[n2];
	strcpy(idvecchunk1,arg[4]);
	strcpy(idvecchunk2,arg[5]);

	// assign atom types for vector computation 
	
	////////////////////////////////////////////////initialize input from run script ////////////////////
	type3 = utils::inumeric(FLERR,arg[3],false, lmp);

	//user defined cutoff for neighbours!
	cutoff_user = utils::numeric(FLERR,arg[6],false, lmp);
	totalsfg = utils::numeric(FLERR,arg[7],false, lmp);

	init();

	// chunk-based data

	nchunk = 1;
	maxchunk = 0;
	allocate();

	

	firstflag = 1;

	//DK: G2, G3 list create
	memory->create(g2g3,atom->natoms,"symf:g2g3");
//if you want to print out some outputs, check g2g3, test1,test2,test3 as an example.


        ostringstream ss;
	ss << me;
	string outfile;
	outfile = "g2g3func-" + ss.str() + ".dat";

        of_g2g3func.open(outfile.c_str(), ofstream::out | ofstream::trunc);

	memory->create(RskappaLst,totalsfg,"g2g3_gen:RskappaLst");
	memory->create(etaLst,totalsfg,"g2g3_gen:etaLst");
	memory->create(sfgnum,6,"g2g3_gen:sfgnum");
	read_input(RskappaLst, etaLst, sfgnum);
}

// ---------------------------------------------------------------------- 
// destructor
// ---------------------------------------------------------------------- 
ComputeG2G3gen::~ComputeG2G3gen()
{
	delete [] idvecchunk1;
	delete [] idvecchunk2;

	memory->destroy(hist);
	memory->destroy(histall);
	memory->destroy(array);
        memory->destroy(g2g3);
	memory->destroy(RskappaLst);
	memory->destroy(etaLst);
}

// ---------------------------------------------------------------------- 
//initialization before a run (optional) (here also called in constructor
// ---------------------------------------------------------------------- 
void ComputeG2G3gen::init()
{
	int icompute1 = modify->find_compute(idvecchunk1);
        int icompute2 = modify->find_compute(idvecchunk2);


	if (icompute1 < 0 || icompute2 < 0)
		error->all(FLERR,"Vector/chunk compute does not exist for compute histo/vecchunk");
	cvecchunk1 = (ComputeVecChunk *) modify->compute[icompute1];
        cvecchunk2 = (ComputeVecChunk *) modify->compute[icompute2];

	if (strcmp(cvecchunk1->style,"vector/chunk") != 0 || strcmp(cvecchunk2->style,"vector/chunk") != 0)
		error->all(FLERR,"Compute histo/vecchunk does not use vector/chunk compute");

	//determine the cutoff for the neighbour list
	double skin = neighbor->skin;
	//mycutneigh = cutoff_user + skin;
	mycutneigh = cutoff_user;

	double cutghost;            // as computed by Neighbor and Comm
	if (force->pair)
		cutghost = MAX(force->pair->cutforce+skin,comm->cutghostuser);
	else
		cutghost = comm->cutghostuser;
	if (mycutneigh > cutghost)
		error->all(FLERR,"Compute histo/vecchunk cutoff exceeds ghost atom range - "
				"use comm_modify cutoff command");

  	neighbor->add_request(this, NeighConst::REQ_FULL);
  	

	// set 1st column of array to bin coords
	for (int i =0 ;i<nbin ;i++){
		array[i][0] = double(i + LOW)/RES;
		array[i][1] = 0.0;
	}
}


// ---------------------------------------------------------------------- 
// init list needed for communciation
// ---------------------------------------------------------------------- 
void ComputeG2G3gen::init_list(int id, NeighList *ptr)
{
  list = ptr;
}



// ---------------------------------------------------------------------- 
// not entirely sure when the setup function is called, re-check
// ---------------------------------------------------------------------- 
void ComputeG2G3gen::setup()
{
	// one-time calculation of per-chunk mass
	// done in setup, so that ComputeChunkAtom::setup() is already called

	//if (firstflag && cvecchunk->idsflag == ONCE) {
	//	compute_array();
	//	firstflag = 0;
	//}
}

// ---------------------------------------------------------------------- 
// here the actual computation is done, here it's array, depending on what's in constructor
// ---------------------------------------------------------------------- 
void ComputeG2G3gen::compute_array()
{
	int indi,indj;
	double massone;
	int i,j,k,m,ii,jj,inum,jnum,ibin,ihisto;
	double histcount,histcountall;  	//count histogram entries  
	int *ilist,*jlist,*numneigh,**firstneigh;

	double cosvec1, cosvec2;

	invoked_array = update->ntimestep;

	// compute chunk/atom assigns atoms to chunk IDs
	// extract ichunk index vector from compute
	// ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms
////////////////////////////////////////////////////////////////////////////////////////////Fix here cvecchunk
	cvecchunk1->compute_array();
	cvecchunk2->compute_array();
	nchunk = cvecchunk1->pu_nchunk;       //check both with cvecchunk1 and 2
	ichunk = cvecchunk1->pu_ichunk;       //check both with cvecchunk1 and 2

	//communicate chunk ids to ghost atoms
	comm->forward_comm(this);


	// compute vector between two atoms of type1 and type2 in chunk

	double **x = atom->x;  		// coordinate vector of local atoms
	int *mask = atom->mask;		// mask if atoms belong to certain group
	int *type = atom->type;		// atom type
	int nlocal = atom->nlocal;  	//total number of local atoms
	tagint *tag = atom->tag;	//global atom id

	//assign molvecall from compute
	double **new_molvec1 = cvecchunk1->pu_molvecall;
	double **new_molvec2 = cvecchunk2->pu_molvecall;
	double rmin, rmax,coss, kappa;
        rmin = cutoff_user-rdiff;
	rmax = cutoff_user;
	double TEMPcf;
	int step=update->ntimestep;
	double xtmp,ytmp,ztmp,delx,dely,delz,r;
	//stuff for the neighbourlists
	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

        for(int i = 0; i < atom->natoms; i++){
		g2g3[i] = 0.0;
	}
	//JR adding some variables for printing

	// zero histograms
	for(i=0;i<1;i++){
		for(j=0;j<nbin;j++){
			hist[i][j] = 0;
		}
	}


	histcount = 0.0;		//initialize histogram counter
	histcountall = 0.0;
	
	int kk = 0;
	for (ii = 0; ii < inum; ii++){		//loop over all atoms on this proc with neighbourlist
		i = ilist[ii];			//index of that atom on this proc
		int id = atom->tag[i];
		if (!(mask[i] & groupbit)){
			continue;  	// only for atoms in the specified group
		}
		indi = ichunk[i]-1;
		if (indi < 0){
			continue;		// only for atoms in chunk
		}
		if(!(type[i] == type3)){                    //change this to type 1 from 3. I dont think I will get different number, but debugging. 
			continue; 	// only for atoms of type 1
		}
		jlist = firstneigh[i];			//neighbours of i
		jnum = numneigh[i]; 			// # neighbours of i
		xtmp = x[i][0];
		ytmp = x[i][1];
		ztmp = x[i][2];
		for (kk = 0; kk < totalsfg; kk++){
			g2g3[kk] = 0.0;
		}
		for(jj = 0; jj < jnum; jj++){	//loop over ii's neighbor
			int flag = 0;
			j = jlist[jj];
			j &= NEIGHMASK;
			if (!(mask[j] & groupbit)) continue;  	// only for atoms in the specified group
			indj = ichunk[j]-1;		
			if (indj < 0) continue;		// only for atoms in chunk
			if(!(type[j] == type3)) continue;	// only atoms of type 1
			delx = xtmp - x[j][0];
			dely = ytmp - x[j][1];
			delz = ztmp - x[j][2];
			r = sqrt(delx*delx + dely*dely + delz*delz);	//finding the distance 
// now indi and indj are the chunk indices of the molvector       (x1*x2 + y1*y2 + z1*z2)/sqrt(x1**2+y1**2+z1**2)/sqrt(x2**2+y2**2+z2**2)   so bekutoru no naiseki tukatte cos motometeru
//change this to CO and NN
			cosvec1 = dot3(new_molvec1[indi],new_molvec1[indj])/(len3(new_molvec1[indi])*len3(new_molvec1[indj]));		//vector distance
			cosvec2 = dot3(new_molvec2[indi],new_molvec2[indj])/(len3(new_molvec2[indi])*len3(new_molvec2[indj]));
			ibin = static_cast<int> (round(cosvec1*RES - LOW));
			if (ibin < 0 || ibin >= nbin) continue;
			hist[0][ibin] += 1.0;
			histcount += 1.0;
			int counter = 0;
				for (int i = 0; i <sfgnum[0];i++)//G2 point
				{
                                        if (abs(r) <= rmin)
						g2g3[counter] += exp(-etaLst[counter]*pow(abs(r)-RskappaLst[counter],2));
					else if ( (rmin < abs(r) && abs(r) <= rmax)){
						TEMPcf = 0.5*(cos((abs(r) - rmin)/(rmax-rmin)*M_PI)+1);
						g2g3[counter] += exp(-etaLst[counter]*pow(abs(r)-RskappaLst[counter],2)) * TEMPcf;
					}
					counter++;
				}

				for (int i = 0; i <sfgnum[1];i++)//G2 point vec 1
				{
	                		if (abs(r) <= rmin)
						g2g3[counter] += exp(-etaLst[counter]*pow(cosvec1-RskappaLst[counter],2));
					else if ( (rmin < abs(r) && abs(r) <= rmax)){
						TEMPcf = 0.5*(cos((r - rmin)/(rmax-rmin)*M_PI)+1);
						g2g3[counter] += exp(-etaLst[counter]*pow(cosvec1-RskappaLst[counter],2)) * TEMPcf;	
					}
					counter++;
				}

				for (int i = 0; i <sfgnum[2];i++)//G2 point vec 2
				{
					if (abs(r) <= rmin)
						g2g3[counter] += exp(-etaLst[counter]*pow(cosvec2-RskappaLst[counter],2));
					else if ( (rmin < abs(r) && abs(r) <= rmax)){
						TEMPcf = 0.5*(cos((r - rmin)/(rmax-rmin)*M_PI)+1);
					       	g2g3[counter] += exp(-etaLst[counter]*pow(cosvec2-RskappaLst[counter],2)) * TEMPcf;
					}
					counter++;
				}

				for (int i = 0; i <sfgnum[3];i++)//G3 point
				{
	                		if (abs(r) <= rmin)
						g2g3[counter] += cos(RskappaLst[counter]*abs(r));
					else if ( (rmin < abs(r) && abs(r) <= rmax)){
						TEMPcf = 0.5*(cos((abs(r) - rmin)/(rmax-rmin)*M_PI)+1);
						g2g3[counter] += cos(RskappaLst[counter]*abs(r))*TEMPcf;
					}
					counter++;
				}	
	
				for (int i = 0; i <sfgnum[4];i++)//G3 point vec 1
				{
					if (abs(r) <= rmin)
						g2g3[counter] += cos(RskappaLst[counter]*(cosvec1));
					else if ( (rmin < abs(r) && abs(r) <= rmax)){
						TEMPcf = 0.5*(cos((r - rmin)/(rmax-rmin)*M_PI)+1);
						g2g3[counter] += cos(RskappaLst[counter]*(cosvec1)) * TEMPcf;
					}
					counter++;
				}
				for (int i = 0; i <sfgnum[5];i++)//G3 point vec 2
				{
					if (abs(r) <= rmin)
						g2g3[counter] += cos(RskappaLst[counter]*(cosvec2));
					else if ( (rmin < abs(r) && abs(r) <= rmax)){
						TEMPcf = 0.5*(cos((r - rmin)/(rmax-rmin)*M_PI)+1);
						g2g3[counter] += cos(RskappaLst[counter]*(cosvec2)) * TEMPcf;
					}
					counter++;
				}
		}		//end of jj neighbors
		for (int k = 0; k < totalsfg; k++){		//print out all symmetry function to g2g3func file. 
			of_g2g3func << " " <<g2g3[k] << "   ";
		}
		of_g2g3func << endl;
	}              //end of ii loop


//////////////////////////////////////////////////////////////////////////////////////////////////
	// sum histograms across procs
	MPI_Allreduce(hist[0],histall[0],1*nbin,MPI_DOUBLE,MPI_SUM,world);
	MPI_Allreduce(&histcount,&histcountall,1,MPI_DOUBLE,MPI_SUM,world);

	// put this in the array
	for(ibin = 0; ibin <nbin; ibin++){
		array[ibin][1] = histall[0][ibin]/histcountall * RES;
	}

			

}


// ---------------------------------------------------------------------- 
//  which information is communicated, packing
// ----------------------------------------------------------------------
int ComputeG2G3gen::pack_forward_comm(int n, int *list, double *buf,
                                          int /*pbc_flag*/, int * /*pbc*/)
{

	int i,j,m;
	
	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		buf[m++] = ichunk[j];
	}
	return m;
}

// ---------------------------------------------------------------------- 
//  which information is communicated, packing
// ----------------------------------------------------------------------
void ComputeG2G3gen::unpack_forward_comm(int n, int first, double *buf)
{
	int i,m,last;
	
	m = 0;
	last = first + n;
	
	for (i = first; i < last; i++) ichunk[i] = buf[m++];
}




/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods insure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */


//----------
//  JR: not sure any of the lock things are needed here...
//----------

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

//void ComputeHistoVecChunk::lock_enable()
//{
//  cvecchunk->lockcount++;
//}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

//void ComputeHistoVecChunk::lock_disable()
//{
//  int icompute = modify->find_compute(idvecchunk);
//  if (icompute >= 0) {
//    cvecchunk = (ComputeVecChunk *) modify->compute[icompute];
//    cvecchunk->lockcount--;
//  }
//}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

//int ComputeHistoVecChunk::lock_length()
//{
//	//nchunk = cchunk->setup_chunks();
//	cvecchunk->compute_array();
//	nchunk = cvecchunk->pu_nchunk;
//	return nchunk;
//}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

//void ComputeHistoVecChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
//{
//  cvecchunk->lock(fixptr,startstep,stopstep);
//}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

//void ComputeHistoVecChunk::unlock(Fix *fixptr)
//{
//  cvecchunk->unlock(fixptr);
//}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeG2G3gen::allocate()
{
	

  maxchunk = nchunk;
  
  
  //assing movecall to array (-> this is accessible via compute)
  //array = molvecall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeG2G3gen::memory_usage()
{
  double bytes = (bigint) maxchunk * 2*3 * sizeof(double);
  bytes += (bigint) maxchunk * 4*3 * sizeof(double);
  return bytes;
}



void ComputeG2G3gen::read_input(double *RskappaLst, double *etaLst,int *sfgnum)
{
	ifstream if_input;	

	string string1;
	getline(if_input,string1);
	
	if_input.open("INPUT.dat",ifstream::in);		//read INPUT.dat
	if(if_input.is_open()){					//checking if the file exist
		do 						//check all the line
		{
			istringstream line(string1);		
			
			if (string1.substr(0,16) =="G2 point        ")	//if the line start with this, take the number of G2 point function to sfgnum
			{
				string test;
				line >> test;			//separate by space
				line >> test;			
				
				line >> test;
				sfgnum[0] = stoi(test);		//store the number of G2 point symmetry function
			}
			if (string1.substr(0,16) =="G2p Rs          ")	//take G2 point Rs parameters
			{
				string test;
				line >> test;
				line >> test;

				for (int i = 0; i<sfgnum[0]; i++)	//loop over number of G2 point function 
				{
					line >> test;
					RskappaLst[i] = stod(test);
				}
			}
			if (string1.substr(0,16) =="G2p eta         ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = 0; i<sfgnum[0]; i++)
				{
					line >> test;
					etaLst[i] = stod(test);		
				}
			}

			if (string1.substr(0,16) =="G2 point vec1   ")	//This is for the first point vector representation you assigned 
			{
				string test;
				line >> test;
				line >> test;
				line >> test;

				line >> test;
				sfgnum[1] = stoi(test);
			}
			if (string1.substr(0,16) =="G2v1 Rs         ")
			{
				string test;
				line >> test;
				line >> test;
				
				for (int i = sfgnum[0]; i<sfgnum[0] + sfgnum[1]; i++)
				{
					line >> test;
					RskappaLst[i] = stod(test);			
				}
			}
			if (string1.substr(0,16) =="G2v1 eta        ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0]; i<sfgnum[0] + sfgnum[1]; i++)
				{
					line >> test;
					etaLst[i] = stod(test);			
				}
			}

			if (string1.substr(0,16) =="G2 point vec2   ")		//This is for the second point vector representation you assigned
			{
				string test;
				line >> test;
				line >> test;
				line >> test;

				line >> test;
				sfgnum[2] = stoi(test);
			}
			if (string1.substr(0,16) =="G2v2 Rs         ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1]; i<sfgnum[0] + sfgnum[1] + sfgnum[2]; i++)
				{	
					line >> test;
					RskappaLst[i] = stod(test);		
				}
			}
			if (string1.substr(0,16) =="G2v2 eta        ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1]; i<sfgnum[0] + sfgnum[1] + sfgnum[2]; i++)
				{
					line >> test;
					etaLst[i] = stod(test);			
				}
			}

			if (string1.substr(0,16) =="G3 point        ")	
			{
				string test;
				line >> test;
				line >> test;
				
				line >> test;
				sfgnum[3] = stoi(test);
			}
			if (string1.substr(0,16) =="G3p kappa       ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1] + sfgnum[2]; i<sfgnum[0] + sfgnum[1] + sfgnum[2] + sfgnum[3]; i++)
				{
					line >> test;
					RskappaLst[i] = stod(test);			
					etaLst[i] = -1;			//G3 function does not need eta. This is why eta is assigned as -1/ 
				}
			}

			if (string1.substr(0,16) =="G3 point vec1   ")
			{
				string test;
				line >> test;
				line >> test;
				line >> test;
				
				line >> test;
				sfgnum[4] = stoi(test);
			}
			if (string1.substr(0,16) =="G3v1 kappa      ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1] + sfgnum[2] + sfgnum[3]; i<sfgnum[0] + sfgnum[1] + sfgnum[2] + sfgnum[3] + sfgnum[4]; i++)
				{
					line >> test;
					RskappaLst[i] = stod(test);		
					etaLst[i] = -1;			
				}
			}

			
			if (string1.substr(0,16) =="G3 point vec2   ")
			{
				string test;
				line >> test;
				line >> test;
				line >> test;
				
				line >> test;
				sfgnum[5] = stoi(test);
			}
			if (string1.substr(0,16) =="G3v2 kappa      ")
			{
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1] + sfgnum[2] + sfgnum[3] + sfgnum[4]; i<sfgnum[0] + sfgnum[1] + sfgnum[2] + sfgnum[3] + sfgnum[4] + sfgnum[5]; i++)
				{	
					line >> test;
					RskappaLst[i] = stod(test);
					etaLst[i] = -1;			
				}
			}
			if (string1.substr(0,16) =="rdiff           ")
			{
				string test;
				line >> test;
				line >> test;
				rdiff = stod(test);
			
			}
                }while (getline(if_input,string1));  //do until end of the file
	}
	else{						//if the file cannot be opened, return the error message. 
		cerr<<"\n\n\t!!ERROR!! Cannot open file "<<endl;
		cerr<<"\tExiting programme..."<<endl<<endl;;
		exit(1);
	}
}
