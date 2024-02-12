#ifndef NNCLASS_H
#define NNCLASS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>


#ifdef MAC_OBLAS
#include <cblas_openblas.h>
#elif OBLAS
#include <cblas.h>
#else
#include <mkl.h>
#endif

#include <omp.h> 
#include <unistd.h>

using namespace std;


// defining some global constants
const int MAX_NODES=51;
const double goodrange[2]={1, 1};  //feature scaling targets (absolute value of range of input and output has to be equal)
double const P_infinity=DBL_MAX*DBL_MAX;
double const N_infinity=-DBL_MAX*DBL_MAX;




// class for all NN variables
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

		double *best;    //best[layers+1][MAX][MAX];  //best weights in the synapses
		double *err;     //err[layers+1][MAX][MAX];  //error gradient
		double *err_priv;
		double *delta;      //delta[layers+2][MAX];  //error in node (delta[0][] should all be empty)
		double *Gdelta;     //Gdelta[layers+2][MAX][MAX];  //dE/dy in node (Gdelta[0][][] should all be empty)
		double *Hv;         //Hv[layers+1][MAX][MAX];  //hessian gradient in error gradient (Herr[0][] should all be empty)
			

		int seed;		//random number seed to initialize weights
		int activation;  // which activation function to use; 0=sigmoid, 1=softmax; also changes loss function accordingly
					//
		double *ZERO;          //reset matrix

};

class SIMvariables{
	private:

	public:
		SIMvariables();
		virtual ~SIMvariables();

		int option;    		// simulation option 1=with new weights, -1=with old weight, -2=with old weights+noise; checking job otherwise
		int niter;		// no. of interations		
		int nprev_iter;  	//
		int nsave;		// how often to save the results
		int noptmax;		// maximum number of steps for optimizer
		

		string trainfile;      // name of input file training data
		string testfile;      // name of input file test data
		string testoutf;
		string testerrf;
		
		string oldweights;     // name of old weights file
		string costfile;       // name of file to save values of cost function
		string newweights;     // name of file for new weights

		int totalpoints;    	//total number of training/test data
		double *dataset;    // contains training/test data
					//
		int nthreads;          // number of ompthreads

};



//-----------------------------------------------------------------------------------------------
//                       function declaration
//-----------------------------------------------------------------------------------------------


void read_input(const string &, NNvariables &, SIMvariables &);     	//reading input file
void train(NNvariables &, SIMvariables &);        				//train the NN
void test(NNvariables &, SIMvariables &);        				//test the NN
double fmincg(double *, int , NNvariables &, SIMvariables &);
void hessproduct(double *, NNvariables &, SIMvariables &);
void grad(double *, NNvariables &, SIMvariables &);
void network(int , double *, NNvariables &);
void hidden(int , int, NNvariables &);
void backprop(int , int, double *, NNvariables &);		// backpropagation
void linemin(double *, double *, NNvariables &, SIMvariables &);               	// line minimisation
double truefunc(double *, NNvariables &, SIMvariables &);
double func(double *, NNvariables &, SIMvariables &);
double h_class(double);  					//activation function
double hprime_class(double);  					//activation function
double softmax(double, double);           //softmax activation function
double errorfunction(double, double);
double errorfunctionprime(double, double);
double crossentropy(double, double, double);        // cross-entropy loss
double crossentropyprime(double, double, double);   // derivative of cross-entropy loss
double accuracy(double, double, double, double,int);
void init(NNvariables &, SIMvariables &);        				//initialise weights
void read_weights(NNvariables &, SIMvariables &);        			//read in old weights





#endif  // NNCLASS_H
