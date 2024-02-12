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

ComputeStyle(G2G3gen/vecchunk,ComputeG2G3gen)

#else

#ifndef LMP_COMPUTE_G2G3GEN_H
#define LMP_COMPUTE_G2G3GEN_H

#include "fix.h"
#include "compute.h"
#include <fstream>
#include <sstream>
using namespace std;
namespace LAMMPS_NS {


	class ComputeG2G3gen : public Compute {
		public:
			char *idvecchunk1;              // fields accessed by other classes
			char *idvecchunk2;
			ComputeG2G3gen(class LAMMPS *, int, char **);
			~ComputeG2G3gen();
			void init();
			void init_list(int, class NeighList *);
			void setup();
			void compute_array();
			void read_input(double *, double *, int *);			// read inout from file

			//to communicate chunk ids to ghost atoms
			int pack_forward_comm(int, int *, double *, int, int *);
			void unpack_forward_comm(int, int, double *);
			double memory_usage();

		private:
			int nchunk,maxchunk;
			int *ichunk;
			int firstflag;
			class ComputeVecChunk *cvecchunk1, *cvecchunk2, *cvecchunk3;

			void allocate();

			//JR start adding some variables
			int me;			//which proc is this
			int nprocs;			//total number of procs

			//stuff for neighbor lists
			class NeighList *list;

			double cutoff_user;    // user-specified cutoff
			double mycutneigh;     // user-specified cutoff + neighbor skin
			
			
			int totalsfg;
			double *RskappaLst;	
			double *etaLst;
			int *sfgnum;
			double rdiff = 0.2;

			//stuff or recording the histogram
			double **hist;         // histogram bins
  			double **histall;      // summed histogram bins across all procs
			int nbin;
			double RES;
			int LOW;

			//DK added for G2 func and G3 func. ofstream is used when we print the data. 
			double *g2g3;
			ofstream of_g2g3func;
			ofstream of_test1func;
			ofstream of_test2func;
			ofstream of_test3func;
			
			// DK addedfrom symf not sure I need all of them. 
			int cutflag;        //may not nessesary
			void init_norm();   //may not nessesary 
			bigint natoms_old;  //may not necessary 
			int type3, type2;
			int nrow, ncol;     //may not necessary
	};

}

#endif
#endif

