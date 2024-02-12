#include "header.h"

#ifndef READINPUT_H
#define READINPUT_H

using namespace std;

// read in number of particles lammps
void get_natoms_lammps(const string&, int &);
// read in particle file lammps
void readParticleFile_lammps(const string&, double **, int&, double *);
void readParticleFile_lammps_aq6(const string&, double **, int&, double *, double *);
void readParticleFile_lammps_allq(const string&, double **, int&, double *, double **);
void readParticleFile_lammps_vec(const string&, double **,  int&, int *, int *, double *);
void readParticleFile_lammps_vec_tri(const string&, double **,  int&, int *, int *, double *);
double minx(double, double);
double maxx(double, double);
double miny(double);
double maxy(double);

#endif // READINPUT_H
