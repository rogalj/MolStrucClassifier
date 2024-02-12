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

ComputeStyle(nnoutgen/atom,ComputeNNoutgen)

#else

#ifndef LMP_COMPUTE_NNOUT_gen_H
#define LMP_COMPUTE_NNOUT_gen_H
#include "compute.h"
#include <fstream>
#include <sstream>
#include "../Dafed.h"
using namespace std;

namespace LAMMPS_NS {

class ComputeNNoutgen : public Compute {
 public:
  ComputeNNoutgen(class LAMMPS *, int, char **);
  ~ComputeNNoutgen();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
//  void set_arrays(int);
//  void refresh();
  double memory_usage();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);        
  void unpack_reverse_comm(int, int *, double *);

 private:
  int nmax;
  double **NNout;
//  char *id_fix;
//  class FixStore *fix;

//  int refreshflag,ivar,nvmax;    // refresh option is enabled
    int nvmax;
//  char *rvar;                    // for incremental dumps
  double *varatom;

  int cutflag;
  double cutoff_user;    // user-specified cutoff
  double mycutneigh;
  class NeighList *list; // half neighbor list
  bigint natoms_old;
  int type1;             //atom type to be used in computing the symmetry function
  int nrow, ncol;
  int me;                        //which proc is this
  int nprocs;
  
  
  int totalsfg;
  double *RskappaLst;	
  double *etaLst;
  int *sfgnum;
  double rdiff = 0.2;

  DAFED::CParameter parameter;
  DAFED::NNvariables NNvar;
  DAFED::CNeighvariables lmpneigh;
  DAFED::CAtom *molecules;                //data structure used within the library        
  DAFED::Cmolpoint *molec;
  
  void read_input(double *, double *,int *,DAFED::CParameter &);
  int nvalues,peratom_freq;
  void get_atominfo(int);
//  double Qall[6];
//  double Qlocal[6];
//  DAFED::EXvariables exvar[2];
  double Rskappa[24]={6.16,6.28,6.76,6.88,0.36,0.08,0.36,0.28,-0.64,-0.36,0.88,1.0,2.5,4.54,4.9,6.22,2.5,3.58,4.78,8.26,2.50,8.12,8.24,8.36};
  double eta[24]={2.44,2.68,1.0,1.0,1.0,1.0,1.12,6.76,3.28,3.28,3.28,3.28,1,1,1,1,1,1,1,1,1,1,1,1};
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for compute displace/atom does not exist

UNDOCUMENTED

E: Compute displace/atom variable is not atom-style variable

UNDOCUMENTED

E: Could not find compute displace/atom fix ID

Self-explanatory.

*/
