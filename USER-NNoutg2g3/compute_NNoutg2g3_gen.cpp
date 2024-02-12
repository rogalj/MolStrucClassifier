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

#include "compute_NNoutg2g3_gen.h"
#include <cmath>
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "modify.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

#include <iostream>
#include <iomanip>
#include "string.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;


/* ---------------------------------------------------------------------- */

ComputeNNoutg2g3gen::ComputeNNoutg2g3gen(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  NNout(nullptr), 
  list(nullptr),
  molecules(),
  molec(),
  nvalues(0),
  lmpneigh(),
  parameter(),
  of_g2g3error(NULL),
  of_g2g3output(NULL),
  of_g2g3good(NULL),
  type1(0)
{
  if (narg !=6) error->all(FLERR,"Illegal compute NNoutgen/atom command");
  if (!atom->tag_enable) error->all(FLERR,"compute symf requires atom tags");
  peratom_flag = 1;
//  size_peratom_cols = 6;
  //create_attribute = 1;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  // optional args
  //type1 = force->inumeric(FLERR,arg[3]);	//carbon atom old code
  type1 = utils::inumeric(FLERR,arg[3],false,lmp);	//carbon atom
  totalsfg = utils::numeric(FLERR,arg[4],false, lmp);
  percent_error = utils::numeric(FLERR,arg[5],false, lmp);
  parameter.nsfg = totalsfg;
//  refreshflag = 0;
//  rvar = nullptr;
  

  cutflag = 1;
  int iarg;
  iarg = 4;
  //DK:checking the cutoff
//  if (strcmp(arg[iarg],"cutoff") != 0) error->all(FLERR,"Illegal compute symf command, require cutoff");
  //cutoff_user = force->numeric(FLERR,arg[iarg+1]);
//  cutoff_user = utils::numeric(FLERR,arg[iarg+1],false,lmp);
//  if (cutoff_user <= 0.0) cutflag = 0;
//  else cutflag = 1;
  
/*  
  iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"refresh") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute NNout/atom command");
      refreshflag = 1;
      delete [] rvar;
      int n = strlen(arg[iarg+1]) + 1;
      rvar = new char[n];
      strcpy(rvar,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal compute NNout/atom command");
  }
*//*
  // error check

  if (refreshflag) {
    ivar = input->variable->find(rvar);
    if (ivar < 0)
      error->all(FLERR,"Variable name for compute NNout/atom does not exist");
    if (input->variable->atomstyle(ivar) == 0)
      error->all(FLERR,"Compute NNout/atom variable "
                 "is not atom-style variable");
  }
*/
  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group
/*
  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[6];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "peratom";
  newarg[4] = (char *) "1";
  newarg[5] = (char *) "3";
  modify->add_fix(6,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;
*/
  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

/*  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;

    double **x = atom->x;
    int *mask = atom->mask;
    imageint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
    	if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
    	else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }
*/
  // per-atom displacement array

  nmax = atom->nmax;
  int natoms=int(atom->natoms);
  DAFED::read_weight(NNvar);
/*  parameter.rmin0 = cutoff_user-0.2;
  parameter.rmax0 = cutoff_user;
  parameter.center = 3;
  parameter.natm = 8;
//  parameter.nsfg=24;
  parameter.nsfg2CO=4;    //no. of G3 symmetry functions
  parameter.nsfg2NN=4;    //no. of G3 symmetry functions
  parameter.nsfg3CO=4;    //no. of G3 symmetry functions
  parameter.nsfg3NN=4;    //no. of G3 symmetry functions
  parameter.nsfg2point=4; //no. of G3 symmetry functions
  parameter.nsfg3point=4; //no. of G3 symmetry functions
  parameter.nmol = natoms/parameter.natm; //no. of total mol in the system
//  parameter.center = 3;
*/  
  parameter.nex = 2;
  molecules = new DAFED::CAtom[nmax];
  molec = new DAFED::Cmolpoint[nmax];
//  parameter.nsfg = 24;
/*  parameter.nnout = 6;
  parameter.COvectype[0] = 3;
  parameter.COvectype[1] = 2;
  parameter.NNvectype[0] = 4;
  parameter.NNvectype[1] = 5;
*/  



  parameter.RskappaLst = new double[parameter.nsfg];
  parameter.etaLst = new double[parameter.nsfg];
  memory->create(RskappaLst,totalsfg,"g2g3_gen:RskappaLst");
  memory->create(etaLst,totalsfg,"g2g3_gen:etaLst");
  memory->create(sfgnum,6,"g2g3_gen:sfgnum");
  read_input(RskappaLst, etaLst, sfgnum,parameter);
  size_peratom_cols = parameter.nnout;
//  read_input(RskappaLst, etaLst, sfgnum);
/*  for (int i = 0; i < totalsfg; i++)
  {
  	cout << RskappaLst[i] << endl;
  	cout << parameter.RskappaLst[i] << endl;
  }
  for (int i = 0; i < totalsfg; i++)
  {
	cout << etaLst[i] << endl;
  	cout << parameter.etaLst[i] << endl;
  }
  for (int i = 0; i < 6; i++)
  {
  	cout << sfgnum[i] << endl;
  }*/
  parameter.nmol = natoms/parameter.natm; //no. of total mol in the system
/*  cout <<"parameter.natm " << parameter.natm << endl;
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
*/


  for(int i=0;i<nmax;i++){
	  molecules[i].sfg = new double[parameter.nsfg];
	  molecules[i].NNout = new double[NNvar.output];
	  molecules[i].NNgrad = new double[NNvar.output*NNvar.input];	
	  molec[i].atom_ID = new int[parameter.natm];
  }


//  cout <<"parameter.nnout checking = " << parameter.nnout << endl;
  //for packing 
  comm_forward = 6+parameter.nnout;
  comm_reverse = parameter.nnout;
  memory->create(NNout,nmax,parameter.nnout,"NNoutgen/atom:NNout");
  array_atom = NNout;
  nmax = nvmax = atom->nmax;
  varatom = nullptr;
  ostringstream ss;
  ss << me;
  string outfileerror;
  string outfilegood;
  string outfileout;
  outfileerror = "g2g3error-" + ss.str() + ".dat";
  outfilegood = "g2g3good-" + ss.str() + ".dat";
  outfileout = "g2g3output-" + ss.str() + ".dat";
  of_g2g3error.open(outfileerror.c_str(), ofstream::out | ofstream::trunc);
  of_g2g3good.open(outfilegood.c_str(), ofstream::out | ofstream::trunc);
  of_g2g3output.open(outfileout.c_str(), ofstream::out | ofstream::trunc);
}

/* ---------------------------------------------------------------------- */

ComputeNNoutg2g3gen::~ComputeNNoutg2g3gen()
{
  // check nfix in case all fixes have already been deleted

  //if (modify->nfix) modify->delete_fix(id_fix);

  //delete [] id_fix;
  memory->destroy(NNout);
  //delete [] rvar;
  memory->destroy(varatom);
  for(int i=0;i<nmax;i++){
	  delete [] molecules[i].sfg;
	  delete [] molecules[i].NNout;
	  delete [] molecules[i].NNgrad;
  	  delete [] molec[i].atom_ID;	  
  }
  delete [] molecules;
  delete [] molec;
}

/* ---------------------------------------------------------------------- */

void ComputeNNoutg2g3gen::init()
{
  // set fix which stores original atom coords
/* 	//DK:commented out
  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute NNoutg2g3/atom fix ID");
  fix = (FixStore *) modify->fix[ifix];

  if (refreshflag) {
    ivar = input->variable->find(rvar);
    if (ivar < 0)
      error->all(FLERR,"Variable name for compute NNout/atom does not exist");
  }
*/
  if (!force->pair)
	  error->all(FLERR,"Compute NNout requires a pair style be defined "
			  " cutoff specified");

  if (cutflag) {
	  double skin = neighbor->skin;	
	  mycutneigh = cutoff_user + skin;

	  double cutghost;            //  as computed by Neighbor and Comm
	  if (force->pair)
		  cutghost = MAX(force->pair->cutforce+skin,comm->cutghostuser);	
	  else
		  cutghost = comm->cutghostuser;
	  if (mycutneigh > cutghost)
		  error->all(FLERR,"Compute symf cutoff exceeds ghost atom range - "
				  "use comm_modify cutoff command");    // JR: user cutoff larger than neighbour list cutoff
	  if (force->pair && mycutneigh < force->pair->cutforce + skin)
		  if (comm->me == 0)
			  error->warning(FLERR,"Compute sym cutoff less than neighbor cutoff - "
					  "forcing a needless neighbor list build");	
	  //JR: this is the bin width
	  //delr = cutoff_user / nbin;	  
  }
  /*for(int i = 0; i < nrow; i++){
	  for(int j = 0; j< ncol; j++){
		  array[i][j] = 0.0;
	  }	
  }*/
//Old Code
//  natoms_old = atom->natoms;
//  dynamic = group->dynamic[igroup];
//  if (dynamic_user) dynamic = 1;
//  int irequest = neighbor->request(this,instance_me);
//  neighbor->requests[irequest]->pair = 0;
//  neighbor->requests[irequest]->compute = 1;
//  neighbor->requests[irequest]->half = 0;
//  neighbor->requests[irequest]->full = 1;
//  neighbor->requests[irequest]->occasional = 1;
//  if (cutflag) {
//	  neighbor->requests[irequest]->cut = 1;
//	  neighbor->requests[irequest]->cutoff = mycutneigh;
//  }

  //auto req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  //if (cutflag) req->set_cutoff(mycutneigh);
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

// ---------------------------------------------------------------------- 
// // init list needed for communciation
// // ---------------------------------------------------------------------- 
void ComputeNNoutg2g3gen::init_list(int /*id*/, NeighList *ptr)
{
	list = ptr;
}
/* ---------------------------------------------------------------------- */

void ComputeNNoutg2g3gen::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  // grow local displacement array if necessary
  int i,j,inum, jnum;
  int *ilist,*jlist,*numneigh,**firstneigh;
  inum = list->inum;              //number of atoms that have a neighbor list
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double box[6];
  tagint *tag = atom->tag;
  tagint *moltag = atom->molecule;
  int natoms=int(atom->natoms);
  int nex = 2;

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
/*  parameter.rmin0 = cutoff_user-0.2;
  parameter.rmax0 = cutoff_user;
  parameter.nsfg = 24;
  parameter.center = 3;
  parameter.natm = 8;
  parameter.nsfg=24;
  parameter.nsfg2CO=4;    //no. of G3 symmetry functions
  parameter.nsfg2NN=4;    //no. of G3 symmetry functions
  parameter.nsfg3CO=4;    //no. of G3 symmetry functions
  parameter.nsfg3NN=4;    //no. of G3 symmetry functions
  parameter.nsfg2point=4; //no. of G3 symmetry functions
  parameter.nsfg3point=4; //no. of G3 symmetry functions
  parameter.nmol = natoms/parameter.natm; //no. of total mol in the system
  parameter.center = 3;
  parameter.nex = 2;
*/  lmpneigh.inum = list->inum;                     //same as local num except using  pair_style hybrid. In this case, number of neigh
  lmpneigh.ilist = list->ilist;                   //get local id for this specific ato
  lmpneigh.numneigh = list->numneigh;             //number of neighbors
  lmpneigh.firstneigh = list->firstneigh;         //list of neighbors for one specific atom
  lmpneigh.nall = atom->nlocal + atom->nghost;    //number of atoms in the processor(includes ghost)
  lmpneigh.me = me;                               //proc number
  lmpneigh.ntot = atom->natoms;
  if (lmpneigh.nall > nmax) {     // re-allocate memory if no. of atoms is larger than nmax!
	  memory->destroy(NNout); 
	  //first delete the variables
	  for(i=0;i<nmax;i++){		//loop over nmax to delete everything
		  delete [] molecules[i].sfg;
		  delete [] molecules[i].NNout;
		  delete [] molecules[i].NNgrad;	
		  delete [] molec[i].atom_ID;
	  }//end loop
	  delete [] molecules;
	  delete [] molec;
	  nmax = lmpneigh.nall;
	  //initialize molecules data structure
	  molecules = new DAFED::CAtom[nmax];
	  molec = new DAFED::Cmolpoint[nmax];
	  //cout <<"parameter.nnout checking = " << parameter.nnout << endl;
	  memory->create(NNout,nmax,parameter.nnout,"NNout/atom:NNout");
	  array_atom = NNout;
	  for(i=0;i<nmax;i++){//loop over new nmax
		  molecules[i].sfg = new double[parameter.nsfg];
		  molecules[i].NNout = new double[NNvar.output];
		  molecules[i].NNgrad = new double[NNvar.output*NNvar.input];
		  molec[i].atom_ID = new int[parameter.natm];
	  }//end loop nmax
  }
  for (i = 0;i<nmax;i++)	//loop over a maximum number of atoms
  {
	for (int outi = 0; outi < parameter.nnout; outi++)	//loop over number of NNoutput
	{					//set up values to be zero for NNout 
		molecules[i].NNout[outi] = 0;	
		NNout[i][outi] = 0;
	}//end of NNout loop
  }//end of nmax loop
  for(i=0;i<lmpneigh.nall;i++){			//loop over owned + ghost atoms
	  for(j=0;j<3;j++){			//loop over xyz 
		  molecules[i].pos[j] = x[i][j];  //get position for each atom    
	  }//end of loop
	  molecules[i].moltype = moltag[i]-1;     //get mol type for each atom. moltag -1 because moltag start from 1 
	  molecules[i].eletype = type[i];         //get element type for each atom
  }//end of nall loop
  DAFED::get_allNeighbourDistances_lmp_new(molecules, molec, lmpneigh, lmpneigh.nall, box, parameter);		//get all neighbor list from here
  //exit(1);
  comm->forward_comm(this);//Need to communicate here to calculate NNgrad
  MPI_Barrier(world);
//  DAFED::get_vec_symmetry_NN_3_lmpneighbour_new(molecules, molec,lmpneigh,box,parameter,Rskappa,eta,NNvar);	//get NNout in here
  DAFED::get_vec_symmetry_NN_3_lmpneighbour_new(molecules, molec,lmpneigh,box,parameter,parameter.RskappaLst,parameter.etaLst,NNvar);
  
  //get molecules.NNout for all the atoms(owned+ghost). You may store in ghost instead of owned, but it will be communicated 
  for (int mi = 0; mi < lmpneigh.lmpmol; mi++)	//loop over lmpmol which is a number of the molecule in the proc
  {
  	for (int atmi = 0; atmi <parameter.natm; atmi++)	//loop over a number of the atoms in each molecule 
	{
		for (int outi = 0; outi < parameter.nnout; outi++){	//loop over number of NNout
			molecules[molec[mi].atom_ID[atmi]].NNout[outi] = molecules[molec[mi].atom_ID[0]].NNout[outi];
		}//end of NNout
  	}//end of atmi loop
  }//end of lmpmol loop

  
  comm->reverse_comm(this);	//all the NNout will be back to owned from ghost  

      


  comm->forward_comm(this);//Need to communicate here to calculate NNgrad
  MPI_Barrier(world);  
//DK added for g2g3 
  double counter[parameter.nnout];	//count which poly is favored
//  comm->forward_comm(this);//Need to communicate here to calculate NNgrad
  bool flag_NNout[lmpneigh.lmpmol];	//check the environment is same as the target NNout. 1 = yes, 0 = no
  for (int mi = 0; mi < lmpneigh.lmpmol; mi++)	//loop over lmpmol which is a number of the molecule in the proc
  {
	for (int n = 0; n < parameter.nnout; n++)	//DK: try to count up the NNout result of neighbor in here. I need to initialize them each time. 
		counter[n] = 0;
  	for (int mj = 0; mj < molec[mi].n_neighbors;mj++)
	{
		int jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index
		int tj =jmol;				//assigning carbon of the neighbor molecule
		int max_index = 0;
		for (int n = 0; n < parameter.nnout; n++)	//DK: try to count up the NNout result of neighbor in here.
		{
			//cout << "molecules[tj].NNout[n]: " << molecules[tj].NNout[n] << endl;
			//if(molecules[molec[mj].atom_ID[0]].NNout[max_index]<molecules[molec[mj].atom_ID[0]].NNout[n])
			if(molecules[tj].NNout[max_index]<molecules[tj].NNout[n])
				max_index=n;
		}	
		counter[max_index] += 1;	//DK: try to get NNout/number_of_layer
	}

	int max_i = 0,max_j = 0;	//max_i ,j for checking highest value for target and neighbor
	for (int n = 0; n < parameter.nnout; n++)	//DK: try to find the index of the highest NNout
	{
		if(molecules[molec[mi].atom_ID[0]].NNout[max_i]<molecules[molec[mi].atom_ID[0]].NNout[n])
			max_i=n;
		cout << "counter[n] " << counter[n] << endl;
		if (counter[max_j] < counter[n])
			max_j=n;
	}
	//DK: in order to get flag to be true, you need to have target and neighbor's poly need to be same AND more than half neighbor is favor to it
/*	cout << "atom_ID: " << molec[mi].atom_ID[0] << endl;
	cout << "counter[max_i]: " << counter[max_i] << "	molec[mi].n_neighbors: " << molec[mi].n_neighbors << endl;
	cout << "counter[max_j]: " << counter[max_j] << endl;
	cout << "counter[max_i]/molec[mi].n_neighbors: " << counter[max_i]/molec[mi].n_neighbors << endl;
*/
	//if (max_i == max_j && counter[max_i]/molec[mi].n_neighbors>0.3)	
	if (counter[max_i]/molec[mi].n_neighbors>percent_error)	
	{
		of_g2g3output << "atom_ID: " << setw(4) << molec[mi].atom_ID[0] << " z_pos " << setw(12) << molecules[molec[mi].atom_ID[0]].pos[2];
		of_g2g3good << "atom_ID: " << setw(4) << molec[mi].atom_ID[0] << " z_pos " << setw(12) << molecules[molec[mi].atom_ID[0]].pos[2];
		for (int outi = 0; outi < parameter.nsfg; outi++)
		{
			of_g2g3output << setw(12) << molecules[molec[mi].atom_ID[0]].sfg[outi];
			of_g2g3good << setw(12) << molecules[molec[mi].atom_ID[0]].sfg[outi];
		}
		for (int outi = 0; outi < parameter.nnout; outi++)
		{
			of_g2g3output << setw(12) << molecules[molec[mi].atom_ID[0]].NNout[outi];
			of_g2g3good << setw(12) << molecules[molec[mi].atom_ID[0]].NNout[outi];
		}
		of_g2g3output << "\n";	
		of_g2g3good << "\n";	

//		cout << "##################GOOD ENV####################" << endl;
	}
	else
	{	
		of_g2g3output << "atom_ID: " << setw(4) << molec[mi].atom_ID[0] << " z_pos " << setw(12) << molecules[molec[mi].atom_ID[0]].pos[2];
		of_g2g3error << "atom_ID: " << setw(4) << molec[mi].atom_ID[0] << " z_pos " << setw(12) << molecules[molec[mi].atom_ID[0]].pos[2];
		for (int outi = 0; outi < parameter.nsfg; outi++)
		{
			of_g2g3output << setw(12) << molecules[molec[mi].atom_ID[0]].sfg[outi];
			of_g2g3error << setw(12) << molecules[molec[mi].atom_ID[0]].sfg[outi];
		}
		for (int outi = 0; outi < parameter.nnout; outi++)
		{
			of_g2g3output << setw(12) << molecules[molec[mi].atom_ID[0]].NNout[outi];
			of_g2g3error << setw(12) << molecules[molec[mi].atom_ID[0]].NNout[outi];
		}
		of_g2g3output << "\n";	
		of_g2g3error << "\n";	
//		cout << "//////////////////WRONG ENV///////////////////" << endl;
	}
  }





  
  //assign NNout from all owned atoms. 
  for (int ti = 0; ti < nlocal; ti++)	//loop over all owned atoms
  {
	for (i = 0; i < parameter.nsfg+parameter.nnout; i++)		//loop over number of NNout
	{
		if (i<parameter.nnout)
			NNout[ti][i] =molecules[ti].NNout[i];
	}//end of NNout loop
  }//end of owned atoms loop
}//end of per_atom
/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   reset per-atom storage values, based on atom-style variable evaluation
   called by dump when dump_modify refresh is set
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */
double ComputeNNoutg2g3gen::memory_usage()
{
  double bytes = (double)nmax* parameter.nnout * sizeof(double);
  bytes += (double)nvmax * sizeof(double);
  return bytes;
}
int ComputeNNoutg2g3gen::pack_forward_comm(int n, int *list, double *buf,int pbc_flag, int *pbc)
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
		for (int nn = 0; nn< parameter.nnout;nn++){
			buf[m++] = molecules[k].NNout[nn];
		}	
	}
	return m;
}
void ComputeNNoutg2g3gen::unpack_forward_comm(int n, int first, double *buf)
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
		for (int nn = 0; nn< parameter.nnout;nn++){
			molecules[i].NNout[nn] += buf[m++];
		}	
	}
}
int ComputeNNoutg2g3gen::pack_reverse_comm(int n, int first, double *buf)
{
	int i,k,num,m,last;
	m = 0;
	last = first + n;
	for (i = first; i < last; i++) {
		for (int k = 0; k< parameter.nnout;k++){
			buf[m++] = molecules[i].NNout[k];
		}	
	}
	return m;	
}
void ComputeNNoutg2g3gen::unpack_reverse_comm(int n, int *list, double *buf)
{
	int i,j,m;
	m = 0;
	for (i = 0; i < n; i++) {
		j = list[i];
		for (int k = 0; k< parameter.nnout;k++){
			molecules[j].NNout[k] += buf[m++];
		}	
	}	
}
void ComputeNNoutg2g3gen::read_input(double *RskappaLst, double *etaLst,int *sfgnum,DAFED::CParameter &parameter)
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
				parameter.nsfg2point = stoi(test);
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
					parameter.RskappaLst[i] = RskappaLst[i];
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
					parameter.etaLst[i] = etaLst[i];
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
				parameter.nsfg2CO = stoi(test);
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
					parameter.RskappaLst[i] = RskappaLst[i];
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
					parameter.etaLst[i] = etaLst[i];
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
				parameter.nsfg2NN = stoi(test);
			}
			if (string1.substr(0,16) =="G2v2 Rs         ")
			{
				//cout << "Checking Rs" << endl;
				string test;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1]; i<sfgnum[0] + sfgnum[1] + sfgnum[2]; i++)
				{	
					line >> test;
					RskappaLst[i] = stod(test);		
					parameter.RskappaLst[i] = RskappaLst[i];
				}
			}
			if (string1.substr(0,16) =="G2v2 eta        ")
			{
				string test;
				//cout << "Checking eta" << endl;
				line >> test;
				line >> test;
				for (int i = sfgnum[0] + sfgnum[1]; i<sfgnum[0] + sfgnum[1] + sfgnum[2]; i++)
				{
				//	cout << "checking para of G2v2 eta";
					line >> test;
					etaLst[i] = stod(test);			
					parameter.etaLst[i] = etaLst[i];
				}
			}

			if (string1.substr(0,16) =="G3 point        ")	
			{
				string test;
				line >> test;
				line >> test;
				
				line >> test;
				sfgnum[3] = stoi(test);
				parameter.nsfg3point = stoi(test);
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
					parameter.RskappaLst[i] = RskappaLst[i];
					etaLst[i] = -1;			//G3 function does not need eta. This is why eta is assigned as -1/ 
					parameter.etaLst[i] = etaLst[i];
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
				parameter.nsfg3CO = stoi(test);
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
					parameter.RskappaLst[i] = RskappaLst[i];
					etaLst[i] = -1;			
					parameter.etaLst[i] = etaLst[i];
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
				parameter.nsfg3NN = stoi(test);
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
					parameter.RskappaLst[i] = RskappaLst[i];
					etaLst[i] = -1;			
					parameter.etaLst[i] = etaLst[i];
				}
			}
			if (string1.substr(0,16) =="cutoff          ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.rmax0 = stod(test);
				
				cutoff_user = stod(test);	
			}
			if (string1.substr(0,16) =="cutoff_diff     ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.rmin0 =  parameter.rmax0 - stod(test);
			}
			if (string1.substr(0,16) =="center_atom     ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.center =  stoi(test);
			}
			if (string1.substr(0,16) =="atm_in_molec    ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.natm =  stoi(test);
			}
			if (string1.substr(0,16) =="num_of_NN_out   ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.nnout =  stoi(test);
			}
			if (string1.substr(0,16) =="eleOfvec1       ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.COvectype[0] =  stoi(test);
				line >> test;
				parameter.COvectype[1] =  stoi(test);
			}
			if (string1.substr(0,16) =="eleOfvec2       ")
			{
				string test;
				line >> test;
				line >> test;
				parameter.NNvectype[0] =  stoi(test);
				line >> test;
				parameter.NNvectype[1] =  stoi(test);
			}
                }while (getline(if_input,string1));  //do until end of the file
	}
	else{						//if the file cannot be opened, return the error message. 
		cerr<<"\n\n\t!!ERROR!! Cannot open file "<<endl;
		cerr<<"\tExiting programme..."<<endl<<endl;;
		exit(1);
	}
}
