/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(vector/chunk,ComputeVecChunk)

#else

#ifndef LMP_COMPUTE_VECTOR_CHUNK_H
#define LMP_COMPUTE_VECTOR_CHUNK_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVecChunk : public Compute {
	public:
	 char *idchunk;              // fields accessed by other classes

	 ComputeVecChunk(class LAMMPS *, int, char **);
	 ~ComputeVecChunk();
	 void init();
	 void setup();
	 void compute_array();

	 void lock_enable();
	 void lock_disable();
	 int lock_length();
	 void lock(class Fix *, bigint, bigint);
	 void unlock(class Fix *);

	 double memory_usage();

	 //JR: some stuff to be accessed from outside for chunks
	 int pu_nchunk;
	 int *pu_ichunk;
	 int lockcount;
	 double **pu_molvecall;
	 int pu_type1;

	private:
	 int nchunk,maxchunk;
	 int firstflag;
	 class ComputeChunkAtom *cchunk;

	 double **molvec,**molvecall;	//vector between 2 atoms (type1->type2) in molecule

	 void allocate();

	 //JR start adding some variables
	 //int nmax;			//max number of owned + ghost atoms on this proc
	 int me;			//which proc is this
	 int nprocs;			//total number of procs
	 double **pos1,**pos1all;  	//position of first atom in vector 
	 double **pos2,**pos2all;  	//position of second atom in vector
	 int type1, type2;		//atom types to mark atoms in molecule for vector

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Chunk/atom compute does not exist for compute com/chunk

Self-explanatory.

E: Compute com/chunk does not use chunk/atom compute

The style of the specified compute is not chunk/atom.

*/
