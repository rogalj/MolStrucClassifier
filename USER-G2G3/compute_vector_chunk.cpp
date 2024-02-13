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

#include "compute_vector_chunk.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "force.h"

//JR: additional header files
#include <iostream>
#include <iomanip>

using namespace LAMMPS_NS;
//JR: add another namespace
using namespace std;

enum{ONCE,NFREQ,EVERY};

// ---------------------------------------------------------------------- 
// constructor
// ---------------------------------------------------------------------- 
ComputeVecChunk::ComputeVecChunk(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  idchunk(NULL), 
  molvec(NULL), molvecall(NULL),
  me(0),
  nprocs(0),
  pos1(0), pos1all(0),
  pos2(0), pos2all(0),
  type1(0), type2(0),
  pu_nchunk(0), pu_ichunk(NULL),
  pu_molvecall(NULL),
  pu_type1(0)
{
	if (narg != 6) error->all(FLERR,"Illegal compute vector/chunk command");

	//This is necessary to properly communicate stuff:
  	if (!atom->tag_enable) error->all(FLERR,"compute vector/chunk requires atom tags");

	//check if this is run on more than one core
	MPI_Comm_rank(world,&me);
  	MPI_Comm_size(world,&nprocs);
	if(nprocs > 1){
		cout<<"Hello, I'm running compute vector/chunk in parallel with nprocs : "<<nprocs<<endl;
		cout<<"This is proc no "<<me
			<<" nlocal = "<<atom->nlocal
			<<endl;
		//error->all(FLERR,"fix dafed only works on SINGLE CORES right now!");
	}


	//in this compute the 'compute_array()' function is called, values storred in variabl 'array'
	array_flag = 1;	
	size_array_cols = 3;
	size_array_rows = 0;
	size_array_rows_variable = 1;
	extarray = 0;

	// ID of compute chunk/atom
	
	int n = strlen(arg[3]) + 1;
	idchunk = new char[n];
	strcpy(idchunk,arg[3]);

	// assign atom types for vector computation 
	// vector type1 -> type2 (=> pos(type2) - pos(type1))
	//type1 = force->inumeric(FLERR,arg[4]);	//old code
	type1 = utils::inumeric(FLERR,arg[4],false, lmp);
	//type2 = force->inumeric(FLERR,arg[5]);	//old code
	type2 = utils::inumeric(FLERR,arg[5],false, lmp);

	init();

	// chunk-based data

	nchunk = 1;
	maxchunk = 0;
	allocate();

	firstflag = 1;
	lockcount = 0;
}

// ---------------------------------------------------------------------- 
// destructor
// ---------------------------------------------------------------------- 
ComputeVecChunk::~ComputeVecChunk()
{
	delete [] idchunk;
	memory->destroy(molvec);
	memory->destroy(molvecall);
	memory->destroy(pos1);
	memory->destroy(pos1all);
	memory->destroy(pos2);
	memory->destroy(pos2all);
}

// ---------------------------------------------------------------------- 
//initialization before a run (optional) (here also called in constructor
// ---------------------------------------------------------------------- 
void ComputeVecChunk::init()
{
	int icompute = modify->find_compute(idchunk);
	if (icompute < 0)
		error->all(FLERR,"Chunk/atom compute does not exist for compute vector/chunk");
	cchunk = (ComputeChunkAtom *) modify->compute[icompute];
	if (strcmp(cchunk->style,"chunk/atom") != 0)
		error->all(FLERR,"Compute com/chunk does not use chunk/atom compute");
}

// ---------------------------------------------------------------------- 
// not entirely sure when the setup function is called, re-check
// ---------------------------------------------------------------------- 
void ComputeVecChunk::setup()
{
	// one-time calculation of per-chunk mass
	// done in setup, so that ComputeChunkAtom::setup() is already called

	if (firstflag && cchunk->idsflag == ONCE) {
		compute_array();
		firstflag = 0;
	}
}

// ---------------------------------------------------------------------- 
// here the actual computation is done, here it's array, depending on what's in constructor
// ---------------------------------------------------------------------- 
void ComputeVecChunk::compute_array()
{
	int index;
	double massone;
	double unwrap[3];

	invoked_array = update->ntimestep;

	// compute chunk/atom assigns atoms to chunk IDs
	// extract ichunk index vector from compute
	// ichunk = 1 to Nchunk for included atoms, 0 for excluded atoms

	nchunk = cchunk->setup_chunks();
	cchunk->compute_ichunk();
	int *ichunk = cchunk->ichunk;

	//things to be accessed public by other classes/functions
	pu_nchunk = nchunk;
	pu_ichunk = ichunk;
	pu_type1 = type1;

	if (nchunk > maxchunk) allocate();
	size_array_rows = nchunk;

	// zero local per-chunk values

	for (int i = 0; i < nchunk; i++){
		molvec[i][0] = molvec[i][1] = molvec[i][2] = 0.0;
		pos1[i][0] = pos1[i][1] = pos1[i][2] = 0.0;
		pos2[i][0] = pos2[i][1] = pos2[i][2] = 0.0;
		molvecall[i][0] = molvecall[i][1] = molvecall[i][2] = 0.0;
	}

	// compute vector between two atoms of type1 and type2 in chunk

	double **x = atom->x;  		// coordinate vector of local atoms
	int *mask = atom->mask;		// mask if atoms belong to certain group
	int *type = atom->type;		// atom type
	int nlocal = atom->nlocal;  	//total number of local atoms

	//JR adding some variables for printing
	int m;
	tagint *tag = atom->tag;
	tagint *moltag = atom->molecule; //this is the tag which molecule the atom belongs to, works with local index!
	int step=update->ntimestep;

	//m = 0;
	//if (me==0){
	//	cout<<"#------------- Step "<<step<<" ------------"<<endl;
	//	cout<<"type1 = "<<type1<<"   type2 = "<<type2<<endl;
	//}


	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit) {  		// only for atoms in the specified group
			index = ichunk[i]-1;		
			if (index < 0) continue;	// only for atoms in chunk
			if ((type[i] != type1) && (type[i] != type2)) continue; 

			// collect positions of atoms of type1 in chunk
			if(type[i] == type1){
				pos1[index][0] = x[i][0];
				pos1[index][1] = x[i][1];
				pos1[index][2] = x[i][2];
			}
			//collect positions of atoms of type2 in chunk
			else if(type[i] == type2){
				pos2[index][0] = x[i][0];
				pos2[index][1] = x[i][1];
				pos2[index][2] = x[i][2];
			}
			else{
  				error->all(FLERR,"Fatal error in compute vector/chunk, requested atom type not found");
			}

			//JR printing out some stuff
			m = tag[i];
			//cout<<"Proc "<<me
			//	//<<"  molecular = "<<molecular
			//	<<"  atom = "<<i
			//	<<"  tag = "<<m
			//	<<"  imol = "<<moltag[i]  // this is the correct index!
			//	//<<"  imol2 = "<<moltag[m]
			//	<<"  chunk = "<<ichunk[i]-1
			//	//<<"  chunk2 = "<<ichunk[m]-1
			//	<<"   type = "<<type[i]
			//	<<"   x = "<<x[i][0]
			//	<<"   y = "<<x[i][1]
			//	<<"   z = "<<x[i][2]
			//	<<endl
			//	<<"   xu = "<<unwrap[0]
			//	<<"   yu = "<<unwrap[1]
			//	<<"   zu = "<<unwrap[2]
			//	<<endl;			

		}
	}

	// exchange information about positions in chunks with other procs
	MPI_Allreduce(&pos1[0][0],&pos1all[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);
	MPI_Allreduce(&pos2[0][0],&pos2all[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);
	
	// pos1all[chunkid][x,y,z] contains positions of atom of type1 in each chunk
	// pos2all[chunkid][x,y,z] contains positions of atom of type2 in each chunk

	//now loop over all chunks to compute CO vector in each molecule
	//this is done by each processor, so quite redundant
	double dx,dy,dz;
	//for (int i = 0; i < nchunk; i++){
	//	dx = pos2all[i][0] - pos1all[i][0];
	//	dy = pos2all[i][1] - pos1all[i][1];
	//	dz = pos2all[i][2] - pos1all[i][2];
	//
	//	domain->minimum_image(dx,dy,dz);
	//
	//	molvecall[i][0] = dx;
	//	molvecall[i][1] = dy;
	//	molvecall[i][2] = dz;
	//}


	//this is a bit better, again looping over atoms on this proc and then selectively computing the vector
	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit) {     		//only for atoms in the specified group
			index = ichunk[i]-1;  
			if (index < 0) continue;   	//only for atoms within the chunk
			if(type[i] == type1){		//only for type1 atoms (so parallel on all procs)
				dx = pos2all[index][0] - pos1all[index][0];
				dy = pos2all[index][1] - pos1all[index][1];
				dz = pos2all[index][2] - pos1all[index][2];

				//cout<<"Chunk "<<index
				//	<<"   Proc "<<me
				//	<<endl;
				//cout<<"BEFORE: dz = "<<dz<<endl;

				domain->minimum_image(dx,dy,dz);

				//cout<<"AFTER: dz = "<<dz<<endl;
				//cout<<endl;
		
				molvec[index][0] = dx;
				molvec[index][1] = dy;
				molvec[index][2] = dz;

			}
		}
	}
	//communicate vector from all procs, molvecall is assigned to the compute array
	MPI_Allreduce(&molvec[0][0],&molvecall[0][0],3*nchunk,MPI_DOUBLE,MPI_SUM,world);
	//assign molvecall to be externally accessible
	pu_molvecall = molvecall;
	//for (int k = 0; k < nchunk; k++){
	//	cout << k << "     " << pu_molvecall[k][0] << "  " <<pu_molvecall[k][1] << pu_molvecall[k][2];
	//}


}

/* ----------------------------------------------------------------------
   lock methods: called by fix ave/time
   these methods insure vector/array size is locked for Nfreq epoch
     by passing lock info along to compute chunk/atom
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   increment lock counter
------------------------------------------------------------------------- */

void ComputeVecChunk::lock_enable()
{
  cchunk->lockcount++;
}

/* ----------------------------------------------------------------------
   decrement lock counter in compute chunk/atom, it if still exists
------------------------------------------------------------------------- */

void ComputeVecChunk::lock_disable()
{
  int icompute = modify->find_compute(idchunk);
  if (icompute >= 0) {
    cchunk = (ComputeChunkAtom *) modify->compute[icompute];
    cchunk->lockcount--;
  }
}

/* ----------------------------------------------------------------------
   calculate and return # of chunks = length of vector/array
------------------------------------------------------------------------- */

int ComputeVecChunk::lock_length()
{
  nchunk = cchunk->setup_chunks();
  return nchunk;
}

/* ----------------------------------------------------------------------
   set the lock from startstep to stopstep
------------------------------------------------------------------------- */

void ComputeVecChunk::lock(Fix *fixptr, bigint startstep, bigint stopstep)
{
  cchunk->lock(fixptr,startstep,stopstep);
}

/* ----------------------------------------------------------------------
   unset the lock
------------------------------------------------------------------------- */

void ComputeVecChunk::unlock(Fix *fixptr)
{
  cchunk->unlock(fixptr);
}

/* ----------------------------------------------------------------------
   free and reallocate per-chunk arrays
------------------------------------------------------------------------- */

void ComputeVecChunk::allocate()
{
  memory->destroy(molvec);
  memory->destroy(molvecall);
	
  memory->destroy(pos1);
  memory->destroy(pos1all);
  memory->destroy(pos2);
  memory->destroy(pos2all);


  maxchunk = nchunk;
  memory->create(molvec,maxchunk,3,"vector/chunk:molvec");
  memory->create(molvecall,maxchunk,3,"vector/chunk:molvecall");
  
  memory->create(pos1,maxchunk,3,"vector/chunk:pos1");
  memory->create(pos1all,maxchunk,3,"vector/chunk:pos1all");
  memory->create(pos2,maxchunk,3,"vector/chunk:pos2");
  memory->create(pos2all,maxchunk,3,"vector/chunk:pos2all");
  
  //assing movecall to array (-> this is accessible via compute)
  array = molvecall;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeVecChunk::memory_usage()
{
  double bytes = (bigint) maxchunk * 2*3 * sizeof(double);
  bytes += (bigint) maxchunk * 4*3 * sizeof(double);
  return bytes;
}
