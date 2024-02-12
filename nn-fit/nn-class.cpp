#include "nn-class.h"


int main(int argc, char *argv[])
{
	NNvariables NNvar;	//all variables for the NN
	SIMvariables SIMvar;	//all simulation variables
	int i;
	string filename="";

	//get input file name from command line input
	if (argc == 2)
	{
		filename = string(argv[1]);
	}
	else{
		cout<<"Programme must be called with 1 arguments: filname!"<<endl;
		cout<<"Exiting programme..."<<endl;
		exit(1);
	}

	//read the input file
	read_input(filename,NNvar,SIMvar);

	// print out some simulation information
	cout<<"------------------------------------------------------------\n";
	cout<<"               Starting NN classification code              \n";
	cout<<"------------------------------------------------------------\n\n";


	cout<<"Random number seed = "<<NNvar.seed<<endl;
	cout<<"Network Architecture:\n";
	cout<<"Input: "<<NNvar.arch[0]<<endl;
	for(i=0;i<NNvar.layers;i++){
		cout<<"Hidden layer "<<i<<": "<<NNvar.arch[i+1]<<endl;
	}
	cout<<"Output: "<<NNvar.arch[NNvar.layers+1]<<endl;
	if(NNvar.activation == 0){
		cout<<"Output activation function: sigmoid"<<endl;
	}
	else if(NNvar.activation == 1){
		cout<<"Output activation function: softmax"<<endl;
	}
	else{
		cout<<"\nError, something is wrong with value of NNvar.activation:" << NNvar.activation << endl;
		cout<<"Exiting programme...\n";
		exit(1);
	}
	cout<<endl;

	cout<<"Scale input data: "<<goodrange[0]<<endl; 
	cout<<"Scale output data: "<<goodrange[1]<<endl;

	cout<<"\nRunning parallel on "<<SIMvar.nthreads<<" cores\n\n";

	cout<<"Option "<<SIMvar.option<<": ";
	switch(SIMvar.option){
		case  1: cout<<"Running training job with new weights\n"; break;
		case -1: cout<<"Running training job with old weights\n"; break;
		//case -2: cout<<"Running training job with old weights + noise\n"; break;
		case 0:  cout<<"Running test job\n"; break;
		default: cout<<"Not doing anything, choose a different run option: -1/1 (training) or 0 (testing)\n";
			 cout<<"Exiting programme...\n";
			 exit(1);
			 break;
	}
	cout<<endl;

	if(SIMvar.option != 0){
		cout<<"Training set file: '"<<SIMvar.trainfile<<"'"<<endl;
		cout<<"Weights output file: '"<<SIMvar.newweights<<"'"<<endl;
		cout<<endl;

		cout<<"Interations: "<<SIMvar.niter<<endl;
		cout<<"Save frequency: "<<SIMvar.nsave<<endl;
		cout<<"Max optimizer steps: "<<SIMvar.noptmax<<endl;
		cout<<"Cost function output file: '"<<SIMvar.costfile<<"'"<<endl;
	}
	else{
		cout<<"Test set file: '"<<SIMvar.testfile<<"'"<<endl;
		cout<<"Weigths file: '"<<SIMvar.oldweights<<"'"<<endl;
		cout<<"Test output: '"<<SIMvar.testoutf<<"'"<<endl;
		cout<<"Test error output: '"<<SIMvar.testerrf<<"'"<<endl;
	}
	cout<<endl;



	
	//initializing arrays 
	NNvar.w        = new double [(NNvar.layers+1)*MAX_NODES*MAX_NODES]; for(i=0; i<(NNvar.layers+1)*MAX_NODES*MAX_NODES; i++)NNvar.w[i]=0;
	NNvar.best     = new double [(NNvar.layers+1)*MAX_NODES*MAX_NODES]; for(i=0; i<(NNvar.layers+1)*MAX_NODES*MAX_NODES; i++)NNvar.best[i]=0;
	NNvar.err      = new double [(NNvar.layers+1)*MAX_NODES*MAX_NODES]; for(i=0; i<(NNvar.layers+1)*MAX_NODES*MAX_NODES; i++)NNvar.err[i]=0;
	NNvar.nodes    = new double [SIMvar.nthreads*(NNvar.layers+2)*MAX_NODES];  for(i=0; i<SIMvar.nthreads*(NNvar.layers+2)*MAX_NODES; i++)NNvar.nodes[i]=0;
	NNvar.delta    = new double [SIMvar.nthreads*(NNvar.layers+2)*MAX_NODES];  for(i=0; i<SIMvar.nthreads*(NNvar.layers+2)*MAX_NODES; i++)NNvar.delta[i]=0;
	NNvar.Gnodes   = new double [SIMvar.nthreads*MAX_NODES*(NNvar.layers+2)*MAX_NODES];  for(i=0; i<SIMvar.nthreads*MAX_NODES*(NNvar.layers+2)*MAX_NODES; i++)NNvar.Gnodes[i]=0;
	NNvar.Gdelta   = new double [SIMvar.nthreads*MAX_NODES*(NNvar.layers+2)*MAX_NODES];  for(i=0; i<SIMvar.nthreads*MAX_NODES*(NNvar.layers+2)*MAX_NODES; i++)NNvar.Gdelta[i]=0;
	NNvar.Hv       = new double [(NNvar.layers+1)*MAX_NODES*MAX_NODES];
	NNvar.err_priv = new double [SIMvar.nthreads*(NNvar.layers+1)*MAX_NODES*MAX_NODES];  for(i=0; i<SIMvar.nthreads*(NNvar.layers+1)*MAX_NODES*MAX_NODES; i++)NNvar.err_priv[i]=0;

	NNvar.stretch  = new double[NNvar.input+NNvar.output];     for(i=0; i<NNvar.input+NNvar.output; i++)NNvar.stretch[i] = 1.0;
	NNvar.translation = new double[NNvar.input+NNvar.output];

	NNvar.ZERO           = new double [(NNvar.layers+1)*MAX_NODES*MAX_NODES]; for(i=0; i<(NNvar.layers+1)*MAX_NODES*MAX_NODES;i++) NNvar.ZERO[i]=0;

	// initialise weights and read training/test data
	init(NNvar,SIMvar);

	// do the training
	if(SIMvar.option != 0){
		train(NNvar,SIMvar);
	}

	// do testing
	if(SIMvar.option == 0){
		test(NNvar,SIMvar);
	}

	return 0;
}


//-----------------------------------------------------------------------------------------------
//                          train the NN
//-----------------------------------------------------------------------------------------------
void train(NNvariables &NNvar, SIMvariables &SIMvar)
{
	double *change;      //change in the weights
	double pre = 1e10;
	double cost = 0.0;
	clock_t t1 = clock();
	int times;
	int i,j,k;

	ofstream f_cost;
	ofstream f_weights;

	f_cost.open(SIMvar.costfile.c_str(), ofstream::out);

	change = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];

	for(times=0; times<SIMvar.niter; times++){
		cost = fmincg(change, times+SIMvar.nprev_iter, NNvar, SIMvar);

		cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1, change, 1, 1, NNvar.best, 1);

		printf("%10d %15e | ETA: %5d seconds||  %15e %15e\n", times+1+SIMvar.nprev_iter, 
			((double)((double)(clock()-t1)/CLOCKS_PER_SEC)/((double)times+1)),
			(int)((double)(SIMvar.niter-times-1)*((double)((double)(clock()-t1)/CLOCKS_PER_SEC)/((double)times+1))),
			cost, pre-cost);
		pre = cost;

		if((times+1+SIMvar.nprev_iter)%5==0){
			// print out cost function
			f_cost << setw(6) << times+1+SIMvar.nprev_iter << "  " << std::scientific<< cost << endl;
		}
		if((times+1+SIMvar.nprev_iter)%SIMvar.nsave == 0){
			//print out weights
			f_weights.open(SIMvar.newweights.c_str(), ofstream::out);
			//architecture
			f_weights << NNvar.layers << endl;
			for(i=0; i<NNvar.layers+2; i++){
				f_weights << NNvar.arch[i] << " ";
			}
			f_weights << endl;
			//activation for output layer
			f_weights << NNvar.activation << endl;
			//translation and stretching
			for(i=0; i<NNvar.input+NNvar.output; i++){
				f_weights << std::scientific << setw(20) << NNvar.translation[i] << " " << setw(20) << NNvar.stretch[i] << " ";
			}
			f_weights << endl;
			// number of training interations
			f_weights << times+1+SIMvar.nprev_iter << endl;

			//weights
			for(i=0; i<NNvar.layers+1; i++){
				for(j=0; j<=NNvar.arch[i]; j++){
					for(k=1; k<=NNvar.arch[i+1]; k++){
						f_weights << std::scientific << setw(15) << NNvar.w[i*MAX_NODES*MAX_NODES+j*MAX_NODES+k] << " ";
					}
					f_weights << endl;
				}
				f_weights << endl << endl;
			}
			f_weights.close();
		}
	}
	f_cost.close();
}

//-----------------------------------------------------------------------------------------------
//                             run NN on test dataset
//-----------------------------------------------------------------------------------------------
void test(NNvariables &NNvar, SIMvariables &SIMvar){
	ofstream f_testout;
	ofstream f_testerr;

	double cost1, cost2, aux;
	int i,n;

	double training[NNvar.input+NNvar.output];  //temporary training set
	double testing[SIMvar.nthreads*NNvar.input];  //given inputs

	double sumexpx;
	double xvalmax;

	cost1 = 0.0;
	cost2 = 0.0;
	sumexpx = 0.0;
	xvalmax = 0.0;

	f_testout.open(SIMvar.testoutf.c_str(), ofstream::out);
	f_testerr.open(SIMvar.testerrf.c_str(), ofstream::out);

	if(f_testout.is_open() && f_testerr.is_open()){
		for(n=0; n<SIMvar.totalpoints; n++){
			memcpy(training, &SIMvar.dataset[n*(NNvar.input+NNvar.output)], sizeof(training));

			for(i=0; i<NNvar.input; i++){
				testing[i] = training[i];
				f_testout << std::scientific << setw(15) << training[i] << " ";
			}
			network(0,testing,NNvar);
			if (NNvar.activation==1){			// collect sum of exponential for softmax function
				xvalmax = -1e10;
			    for (i=1; i<=NNvar.output; i++) xvalmax = max(xvalmax,NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]);	
				sumexpx = 0.0;
				for (i=1; i<=NNvar.output; i++){
					sumexpx += exp(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-xvalmax);
				}
			}
			aux = 0.0;
			for(i=1; i<=NNvar.output; i++){
				aux += accuracy(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i], training[NNvar.input+i-1],xvalmax,sumexpx,NNvar.activation);
			}
			if(aux > 0.0){
				//print to error file
				for(i=0; i<NNvar.input; i++){
					f_testerr << std::scientific << setw(15) << training[i] << " ";
				}
				for(i=1; i<=NNvar.output; i++){
					if (NNvar.activation==0) f_testerr << std::scientific << setw(15) << h_class(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]) << " ";
					else if (NNvar.activation==1) f_testerr << std::scientific << setw(15) << softmax(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-xvalmax,sumexpx) << " ";
					else{
						cout << "\n Something wrong with variable NNvar.activation in test function 1 ... exiting programme \n";
						exit(0);
					}
				}
				for(i=1; i<=NNvar.output; i++){
					f_testerr << std::scientific << setw(15) << training[NNvar.input+i-1] << " ";
				}
				f_testerr << endl;
				cost1 += 1.0;
			}
			for(i=1; i<=NNvar.output; i++){
				if (NNvar.activation==0) f_testout << std::scientific << setw(15) << h_class(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]) << " ";
				else if (NNvar.activation==1) f_testout << std::scientific << setw(15) << softmax(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-xvalmax,sumexpx) << " "; 
				else{
						cout << "\n Something wrong with variable NNvar.activation in test function 2 ... exiting programme \n";
						exit(1);
					}
			}
			for(i=1; i<=NNvar.output; i++){
				f_testout << std::scientific << setw(15) << training[NNvar.input+i-1] << " ";
			}
			f_testout << endl;

			for(i=1; i<=NNvar.output; i++){
				if (NNvar.activation == 0) cost2 += errorfunction(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i], training[NNvar.input+i-1]);
				else if(NNvar.activation == 1) cost2 += crossentropy(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-xvalmax, training[NNvar.input+i-1],sumexpx);
				else{
					cout << "\n Something wrong with variable NNvar.activation in test function 3 ... exiting programme \n";
					exit(1);	
				} 
			}

			if(n%10000 == 0){
				cout << n << " L1 accuracy: "<< std::scientific << setw(15) << cost1/n << " \n";
			}

		}
		cout<<"L1 accuracy: "<< std::fixed << setprecision(3) << (1-cost1/SIMvar.totalpoints)*100 << endl;
		cout<<"Cost function: " << std::scientific << setprecision(-1) << setw(15) << sqrt(cost2/SIMvar.totalpoints) <<" \n";
	}
	else{
		cerr << "\nFatal Error: cannot open test output files '"<<SIMvar.testoutf<<"' and '"<<SIMvar.testerrf<<"' \n";
		cerr << "Exiting programme..."<<endl<<endl;
		exit(1);
	}

	f_testout.close();
	f_testerr.close();


}

//-----------------------------------------------------------------------------------------------
//                      do conjugate gradient step
//-----------------------------------------------------------------------------------------------
double fmincg(double *X, int iter, NNvariables &NNvar, SIMvariables &SIMvar)
{
	double thresh = 1e-30;
	double cost = 0.0;
	double *derivative;
	long double beta = 0.0;
	long double denom = 0.0;
	double pre = 1e10, previous = 1e10;
	double *direction;
	int times;
	double *temp;

	derivative = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];
	direction = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];

	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.ZERO, 0, direction, 1);
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.ZERO, 0, X, 1);


	temp = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];

	for(times=0; times<(NNvar.layers+1)*MAX_NODES*MAX_NODES && times<SIMvar.noptmax; times++){
		// get gradient
		hessproduct(X, NNvar, SIMvar);
		grad(NNvar.best,NNvar,SIMvar);
		cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.err, 1, derivative, 1);
		cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1, NNvar.Hv, 1, 1, derivative, 1);

		//update direction cancelation factor (initial is 0)
		if(times!=0){
			hessproduct(direction, NNvar, SIMvar);
			denom=cblas_ddot((NNvar.layers+1)*MAX_NODES*MAX_NODES, direction, 1, NNvar.Hv, 1);
			if(denom>=-thresh && denom<=thresh)
			{
				delete[] derivative;
				delete[] direction;
				return cost;
			}
			beta=cblas_ddot((NNvar.layers+1)*MAX_NODES*MAX_NODES, derivative, 1, NNvar.Hv, 1);
			beta=beta/denom;
		}
		//update direction
		cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, -1, derivative, 1, beta, direction, 1);

		//line minimize in direction
		linemin(X, direction,NNvar,SIMvar);

		cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, X, 1, temp, 1);
		cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1, NNvar.best, 1, 1, temp, 1);

		cost = func(temp,NNvar,SIMvar);

		printf("%10d | %10d / %10d ||  %15e %15e\n", iter+1, times+1, (NNvar.layers+1)*MAX_NODES*MAX_NODES, cost, pre-cost);

		if((pre-cost)/cost<1e-5 && previous/pre<1e-5){
			delete[] derivative;
			delete[] direction;
			return cost;
		}

		previous = pre-cost;
		pre = cost;
	}

	delete[] derivative;
	delete[] direction;

	return cost;
}

//-----------------------------------------------------------------------------------------------
//                      Hessian product with optimal direction
//-----------------------------------------------------------------------------------------------
void hessproduct(double *direction, NNvariables &NNvar, SIMvariables &SIMvar)
{
	double small = 1e-10;
	double *move;

	move = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.best, 1, move, 1);
	cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, small, direction, 1, 1, move, 1);

	grad(move,NNvar,SIMvar);

	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.err, 1, NNvar.Hv, 1);

	grad(NNvar.best,NNvar,SIMvar);

	cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, -1/small, NNvar.err, 1, 1/small, NNvar.Hv, 1);

	delete[] move;

}

//-----------------------------------------------------------------------------------------------
//                      compute gradient
//-----------------------------------------------------------------------------------------------
void grad(double *weights, NNvariables &NNvar, SIMvariables &SIMvar)
{
	int n, nth=0;
	int i,k;
	
	double training[NNvar.input+NNvar.output];  //temporary training set
	double testing[SIMvar.nthreads*NNvar.input];  //given inputs
	double target[SIMvar.nthreads*NNvar.output];  //given outputs


	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, weights, 1, NNvar.w, 1);
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.ZERO, 0, NNvar.err, 1);
    cblas_dcopy((NNvar.layers+1)*SIMvar.nthreads*MAX_NODES*MAX_NODES, NNvar.ZERO, 0, NNvar.err_priv, 1);

	#pragma omp parallel num_threads(SIMvar.nthreads) //shared(node,err_priv)//lastprivate(err_priv)
	{
	#pragma omp  for private(nth,n,training) schedule(static)
	for(n=0; n<SIMvar.totalpoints; n++){
		nth = omp_get_thread_num();
		memcpy(training, &SIMvar.dataset[n*(NNvar.input+NNvar.output)], sizeof(training));
		for(i=0; i<NNvar.input; i++){
			testing[nth*NNvar.input+i]=training[i];
		}
		for(i=0; i<NNvar.output;i++){
			target[nth*NNvar.output+i]=training[NNvar.input+i];
		}
		network(nth,testing,NNvar);

		backprop(nth,1,target,NNvar);
	}
        #pragma omp barrier
	}

	for(k=0; k<SIMvar.nthreads;k++){
		cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1, &NNvar.err_priv[k*(NNvar.layers+1)*MAX_NODES*MAX_NODES], 1, 1, NNvar.err, 1);
	}
	cblas_dscal((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1.0/SIMvar.totalpoints, NNvar.err, 1);  //average gradient
}

//-----------------------------------------------------------------------------------------------
//                       going through the network
//-----------------------------------------------------------------------------------------------
void network(int nth, double *testing, NNvariables &NNvar)
{
	int i;
	//int lay;

	NNvar.nodes[nth*(NNvar.layers+1)*MAX_NODES+0] = 1.0;
	for(i=1; i<=NNvar.input;i++){
		NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+i] = (testing[nth*NNvar.input+i-1]-NNvar.translation[i-1])*NNvar.stretch[i-1];
	}
	//lay = 1;   // which hidden layer
	hidden(nth,1,NNvar);

	for(i=1; i<=NNvar.output; i++){  //output layer without activation
		NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+(NNvar.layers+1)*MAX_NODES+i] = 
			(NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+(NNvar.layers+1)*MAX_NODES+i])/NNvar.stretch[i-1+NNvar.input]+NNvar.translation[i-1+NNvar.input];
	}


}


//-----------------------------------------------------------------------------------------------
//                 go through hidden layers
//-----------------------------------------------------------------------------------------------
void hidden(int nth, int lay, NNvariables &NNvar)  // this only goes over hidden layers
{
	double *temp;
	int i;

	temp = new double[MAX_NODES];
	cblas_dcopy(MAX_NODES, NNvar.ZERO, 1, temp, 1);
	temp[0] = 1.0;
	
	for(i=1; i<=NNvar.arch[lay-1] && lay==1; i++){
		temp[i] = NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+(lay-1)*MAX_NODES+i];
	}
	for(i=1; i<=NNvar.arch[lay-1] && lay!=1; i++){  // only for hidden layers
		temp[i] = h_class(NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+(lay-1)*MAX_NODES+i]);
	}

	cblas_dgemv(CblasRowMajor, CblasTrans, NNvar.arch[lay-1]+1, NNvar.arch[lay]+1, 1,
		&NNvar.w[(lay-1)*MAX_NODES*MAX_NODES], MAX_NODES, temp, 1, 0, &NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES], 1);

	NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES]=1.0;
	
	delete[] temp;
	if(lay < NNvar.layers+1){
		hidden(nth,lay+1,NNvar);
	}

}

//-----------------------------------------------------------------------------------------------
//                    backpropagation
//-----------------------------------------------------------------------------------------------
void backprop(int nth, int lay, double *target, NNvariables &NNvar)  // this goes over all layers, also the output
{
	int i;
	double sumexpx;   //sum over exponential of outputs for softmax
	double xvalmax;
	double predict;
	sumexpx = 0.0;
	xvalmax = 0.0;
	
	if(lay+1 <= NNvar.layers+1) backprop(nth,lay+1,target,NNvar);

	//get sum of exponentials for softmax
	if (lay == NNvar.layers+1 && NNvar.activation == 1){
		xvalmax = -1e10;
		for(i=1; i<=NNvar.arch[lay]; i++){
			predict = (NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input];
			xvalmax = max(xvalmax,predict);
		}
		sumexpx = 0.0;
		for(i=1; i<=NNvar.arch[lay]; i++){
			predict = (NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input];
			sumexpx += exp(predict-xvalmax);
		}
	}
	if(lay == NNvar.layers+1){  // for the output layer
		for(i=1; i<=NNvar.arch[lay]; i++){
			if (NNvar.activation == 0){
				NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i] = 
					errorfunctionprime((NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input],
							(target[nth*NNvar.output+i-1]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input]);
			}
			else if (NNvar.activation == 1){  // this should be fine with dL/dO_i
				predict = (NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input]; 
				NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i] = 
					crossentropyprime(predict-xvalmax,
							(target[nth*NNvar.output+i-1]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input],sumexpx);	
			}
			else{
				cout<<"\n Something wrong with variable NNvar.activation in function backprop...exiting programme\n";
				exit(1);
			}
		}
		
		NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES] = 0.0;   //bias does not pass error
	}
	else{  // for all other layers
		cblas_dgemv(CblasRowMajor, CblasNoTrans, MAX_NODES, MAX_NODES, 1,
			&NNvar.w[lay*MAX_NODES*MAX_NODES], MAX_NODES, &NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+(lay+1)*MAX_NODES], 
			1, 0, &NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES], 1);

		for(i=1; i<=NNvar.arch[lay]; i++){
			NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i] = 
				NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i]*hprime_class(NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES+i]);
		}
		
		NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+lay*MAX_NODES]=0.0;    //bias does not pass error
	}
	// calculate error gradient
	
	double j;
	for(i=0; i<=NNvar.arch[lay-1];i++){
		if(i==0) j=1.0;
		else if(lay-1==0) j=NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+(lay-1)*MAX_NODES+i];
		else j=h_class(NNvar.nodes[nth*(NNvar.layers+2)*MAX_NODES+(lay-1)*MAX_NODES+i]);  // this should be fine since it is not the output layer (?) 

		cblas_daxpby(MAX_NODES, j, &NNvar.delta[nth*(NNvar.layers+2)*MAX_NODES+(lay)*MAX_NODES], 1, 1, 
				&NNvar.err_priv[nth*(NNvar.layers+1)*MAX_NODES*MAX_NODES + (lay-1)*MAX_NODES*MAX_NODES+i*MAX_NODES], 1);
	}


}

//-----------------------------------------------------------------------------------------------
//                                 perform line minimisation
//-----------------------------------------------------------------------------------------------
void linemin(double *X, double *direction, NNvariables &NNvar,SIMvariables &SIMvar)
{
	double alpha = 1e0;
	double gamma = 0.5;
	double tau = 0.1;
	double *unit;
	double *here;
	double *start;
	double accounted = 0.0;
	long double magnitude = 0.0;

	unit = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];
	here = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];
	start = new double[(NNvar.layers+1)*MAX_NODES*MAX_NODES];

	magnitude = cblas_ddot((NNvar.layers+1)*MAX_NODES*MAX_NODES, direction, 1, direction, 1);
	magnitude = sqrt(magnitude);
	if(magnitude==0){
		delete[] unit;
		delete[] start;
		delete[] here;
		return;
	}
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, direction, 1, unit, 1);
	cblas_dscal((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1/magnitude, unit, 1);
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, X, 1, start, 1);
	cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, 1, NNvar.best, 1, 1, start, 1);
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, start, 1, here, 1);
	cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, alpha, unit, 1, 1, here, 1);
	
	grad(start,NNvar,SIMvar);
	
	accounted = cblas_ddot((NNvar.layers+1)*MAX_NODES*MAX_NODES, unit, 1, NNvar.err, 1);

	double a = truefunc(start,NNvar,SIMvar);
	double b = truefunc(here,NNvar,SIMvar);
	while(a-b < -alpha*gamma*accounted){
		alpha=alpha*tau;
		if(alpha<1e-10){
			delete[] unit;
			delete[] start;
			delete[] here;
			return;
		}
		cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, start, 1, here, 1);
		cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, alpha, unit, 1, 1, here, 1);
		b = truefunc(here,NNvar,SIMvar);
	}
	cblas_daxpby((NNvar.layers+1)*MAX_NODES*MAX_NODES, alpha, unit, 1, 1, X, 1);


}

//-----------------------------------------------------------------------------------------------
//                                    function without feature scaling?
//-----------------------------------------------------------------------------------------------
double truefunc(double *weights, NNvariables &NNvar, SIMvariables &SIMvar)  //without feature scaling
{
	double cost,costing;
	int i,n;
	double training[NNvar.input+NNvar.output];  //temporary training set
	double testing[SIMvar.nthreads*NNvar.input];  //given input 
	double sumexpx;
	double xvalmax;
	double predict; 
	
	sumexpx = 0.0;
	xvalmax = 0.0;
	cost = 0.0;
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, weights, 1, NNvar.w, 1);
	
	for(n=0; n<SIMvar.totalpoints; n++){
		memcpy(training, &SIMvar.dataset[n*(NNvar.input+NNvar.output)], sizeof(training));
		for(i=0; i<NNvar.input; i++) testing[i]=training[i];
		
		network(0,testing,NNvar);

		if(NNvar.activation == 1){
			xvalmax = -1e10;
			for(i=1; i<=NNvar.output; i++){
				predict = (NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input];
				xvalmax = max(xvalmax,predict); 
			}
			sumexpx = 0.0;
			for(i=1; i<=NNvar.output; i++){
				predict = (NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input];
				sumexpx += exp(predict-xvalmax); 
			}
		}
		
		costing = 0.0;
		for(i=1; i<=NNvar.output; i++){
			if (NNvar.activation == 0){
				costing+=errorfunction((NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input],
					(training[NNvar.input+i-1]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input]);
			}
			else if (NNvar.activation == 1){
				predict = (NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input];
				costing+=crossentropy(predict-xvalmax,
					(training[NNvar.input+i-1]-NNvar.translation[i-1+NNvar.input])*NNvar.stretch[i-1+NNvar.input],sumexpx);		
			}
			else{
				cout<<"\n Error, something wrong with variable NNvar.activation in function truefunc...exiting programme\n";
				exit(1);
			}
		}
		cost+=costing/SIMvar.totalpoints;
	}

	if(cost!=cost || cost==P_infinity || cost==N_infinity)cost=DBL_MAX;
	
	return cost;
}

//-----------------------------------------------------------------------------------------------
//                   compute cost function
//-----------------------------------------------------------------------------------------------
double func(double *weights, NNvariables &NNvar, SIMvariables &SIMvar)
{
	double cost,costing;
	int i,n;
	double training[NNvar.input+NNvar.output];  //temporary training set
	double testing[SIMvar.nthreads*NNvar.input];  //given inputs
	double sumexpx;
	double xvalmax;
	double predict; 
	
	sumexpx = 0.0;
	xvalmax = 0.0;

	cost = 0.0;
	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, weights, 1, NNvar.w, 1);
	
	for(n=0; n<SIMvar.totalpoints;n++){
		memcpy(training, &SIMvar.dataset[n*(NNvar.input+NNvar.output)], sizeof(training));
		
		for(i=0; i<NNvar.input; i++) testing[i]=training[i];
		
		network(0,testing,NNvar);

		if(NNvar.activation == 1){
			xvalmax = -1e10;
			for(i=1; i<=NNvar.output; i++){
				predict = NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i];
				xvalmax = max(xvalmax,predict); 
			}
			sumexpx = 0.0;
			for(i=1; i<=NNvar.output; i++){
				predict = NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i];
				sumexpx += exp(predict-xvalmax); 
			}
		}
		
		costing = 0.0;
		for(i=1; i<=NNvar.output; i++){
			if (NNvar.activation == 0) costing += errorfunction(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i], training[NNvar.input+i-1]);
			else if (NNvar.activation == 1) costing += crossentropy(NNvar.nodes[(NNvar.layers+1)*MAX_NODES+i]-xvalmax, training[NNvar.input+i-1],sumexpx); 
			else{
				cout<<"\n Error, something wrong with variable NNvar.activation in function func...exiting programme\n";
				exit(1);
			}
		}
		
		cost+=costing/SIMvar.totalpoints;
	}
	
	if(cost!=cost || cost==P_infinity || cost==N_infinity) cost=DBL_MAX;
	
	return cost;
}

//-----------------------------------------------------------------------------------------------
//                 activation function
//-----------------------------------------------------------------------------------------------
double h_class(double x)  //sigmoid
{
	double temp;
	temp = 1.0/(1.0+exp(-x));
	return temp;
}

double softmax(double x, double sumexpx)   //softmax
{
	double temp;
	temp = exp(x)/sumexpx;
	return temp; 
}

//-----------------------------------------------------------------------------------------------
//                  derivative of activation function
//-----------------------------------------------------------------------------------------------
double hprime_class(double x)
{ // WATCH OUT FOR INFINITY!!!!!
	double temp;
	temp=(double)1/(((double)1+exp(-x))*((double)1+exp(x))); //h_class(x)*((double)1 - h_class(x)) = 1.0/((1+exp(-x))*(1+exp(x)));
	return temp;
}


//-----------------------------------------------------------------------------------------------
//                  loss function
//-----------------------------------------------------------------------------------------------

double errorfunction(double predict, double actual)  //for two classes
{
	if(actual == 0 ) { return -log((double)1-h_class(predict)); } 
	else             { return -log(h_class(predict)); }
}

double crossentropy(double predict, double actual, double sumexpx)
{
	return (-actual*log(softmax(predict,sumexpx)));
}

//-----------------------------------------------------------------------------------------------
//                 derivative of loss function
//-----------------------------------------------------------------------------------------------
double errorfunctionprime(double predict, double actual)
{
	if(actual == 0 ) { return  hprime_class(predict)/((double)1-h_class(predict)); }
	else             { return -hprime_class(predict)/(h_class(predict)); }
}

double crossentropyprime(double predict, double actual, double sumexpx)
{
	return (softmax(predict,sumexpx) - actual);
}

//-----------------------------------------------------------------------------------------------
//                  check if outputs are consistent
//-----------------------------------------------------------------------------------------------
double accuracy(double predict, double actual, double xvalmax, double sumexpx, int activation)
{
	double tmp;
	tmp = 0.0;

	if(activation == 0) tmp = h_class(predict);
	else if (activation == 1) tmp = softmax(predict-xvalmax,sumexpx);
	else{
		cout << "\n Something wrong with variable activation in function accuracy ... exiting programme\n ";
		exit(1);
	}
	if(tmp > 0.5) { return fabs(1.0-actual); }
	else            { return actual ; }  
}


//-----------------------------------------------------------------------------------------------
//                           initialise weights and read training data
//-----------------------------------------------------------------------------------------------
void init(NNvariables &NNvar, SIMvariables &SIMvar)
{
	int i,j,k;
	double scale;

	if(SIMvar.option == 1){  		//randomly initialize weights
		srand(NNvar.seed);
		for(i=0; i<NNvar.layers+1; i++){
			scale = (double) sqrt(6.0/(double)(NNvar.arch[i]+NNvar.arch[i+1]));
			for(j=0; j<=NNvar.arch[i]; j++){
				for(k=1; k<=NNvar.arch[i+1];k++){
					NNvar.w[i*MAX_NODES*MAX_NODES+j*MAX_NODES+k]=0;
					while(NNvar.w[i*MAX_NODES*MAX_NODES+j*MAX_NODES+k]==0 && j!=0){
						NNvar.w[i*MAX_NODES*MAX_NODES+j*MAX_NODES+k]=scale*((double)(rand()%201-100)/100);
					}
				}
			}
		}
	}

	if(SIMvar.option <= 0){
		//read in weights file
		read_weights(NNvar,SIMvar);
	}

	cblas_dcopy((NNvar.layers+1)*MAX_NODES*MAX_NODES, NNvar.w, 1, NNvar.best, 1);

	// read in training/testing data
	ifstream f_data;
	double range[NNvar.input+NNvar.output][2];

	if(SIMvar.option == 0){
		f_data.open(SIMvar.testfile.c_str(), ifstream::in);
		cout<<"\nReading in test data from '"<<SIMvar.testfile<<"' \n";
	}
	else{
		f_data.open(SIMvar.trainfile.c_str(), ifstream::in);
		cout<<"\nReading in training data from '"<<SIMvar.trainfile<<"' \n";
	}
	
	if (f_data.is_open()){

		string line_str;
		char line_in[1000];
		int runi;

		char *stringpointer;
		int readcolumn;
		char *col[256];

		double trainval;

		
		getline(f_data,line_str);

		SIMvar.totalpoints = 0;


		// read in training/test data line by line
		do{
			strcpy(line_in,line_str.c_str());
			stringpointer = line_in;
			readcolumn = 0;

			// run over characters in each line and separate them into columns
			for(runi=0;runi<1000;runi++){
				if(line_in[runi] == 0){    // if you are at the end of the line
					break;
				}
				if(line_in[runi] == ' ' && line_in[runi+1] != ' '){   // get values separated by spaces
					line_in[runi] = 0;
					col[readcolumn] = stringpointer;
					stringpointer = &(line_in[runi+1]);
					readcolumn++;
				}
			}

			if(readcolumn != NNvar.input+NNvar.output){
				cerr << "\n!! ERROR !! Reading training data\n";
				cerr << "No. of columns does not match no. of input + output values\n";
				cerr << "Values in a line must be separated by spaces\n";
				cerr << "Exiting programme..."<<endl<<endl;
				for(i=0;i<readcolumn;i++){
					cout<<i<<" "<<col[i]<<endl;
				}
				exit(1);
			}

			// only if option > 0
			for(i=0; i<NNvar.input+NNvar.output && SIMvar.option>0 ; i++){
				if(SIMvar.totalpoints == 0) NNvar.translation[i] = 0;
				trainval = atof(col[i]);
				NNvar.translation[i] += trainval;

				if(SIMvar.totalpoints == 0 || range[i][0] > trainval) range[i][0] = trainval;
				if(SIMvar.totalpoints == 0 || range[i][1] < trainval) range[i][1] = trainval;
			}
			
			// check proper range of output values for classification
			for(i=NNvar.input; i<NNvar.input+NNvar.output && SIMvar.option>0 ; i++){
				if(atof(col[i]) != 0.0 && atof(col[i]) != 1.0){
					cerr << "\n!! ERROR !! Reading training/test data\n";
					cerr << "For classification, labels should be either 0 or 1!\n";
					cerr << "Exiting programme...\n";
					exit(1);
				}
			}

			SIMvar.totalpoints += 1;

		}while (getline(f_data,line_str));

		// only if option > 0
		for(i=0; i<NNvar.input && SIMvar.option>0 ; i++){ 
			NNvar.translation[i] = NNvar.translation[i]/SIMvar.totalpoints;
			NNvar.stretch[i] = goodrange[0]/((range[i][1]-range[i][0])/2);
		}
		for(i=NNvar.input; i<NNvar.input+NNvar.output && SIMvar.option>0 ; i++){
			NNvar.translation[i] = 0.0;                          // no translation for classification
			NNvar.stretch[i]=1.0/((range[i][1]-range[i][0]));   // for classification only between 0 and 1
		}

		SIMvar.dataset = new double [(SIMvar.totalpoints+1)*(NNvar.input+NNvar.output)];

	}
	else{
		if(SIMvar.option != 0){
			cerr << "\nFatal Error: cannot open training data file '"<<SIMvar.trainfile<<"' \n";
		}
		else{
			cerr << "\nFatal Error: cannot open test data file '"<<SIMvar.testfile<<"' \n";
		}
		cerr << "Exiting programme..."<<endl<<endl;
		exit(1);
	}

	f_data.close();

	//now read in all data to the dataset array
	if(SIMvar.option == 0){
		f_data.open(SIMvar.testfile.c_str(), ifstream::in);
	}
	else{
		f_data.open(SIMvar.trainfile.c_str(), ifstream::in);
	}
	if (f_data.is_open()){
		int n;
		for (n=0; n<SIMvar.totalpoints; n++){
			for(i=0; i<NNvar.input+NNvar.output;i++){
				f_data >> SIMvar.dataset[n*(NNvar.input+NNvar.output) + i];
			}
		}
	}
	else{
		if(SIMvar.option != 0){
			cerr << "\nFatal Error: cannot open training data file '"<<SIMvar.trainfile<<"' , 2nd time\n";
		}
		else{
			cerr << "\nFatal Error: cannot open test data file '"<<SIMvar.testfile<<"' , 2nd time\n";
		}
		cerr << "Exiting programme..."<<endl<<endl;
		exit(1);
	}
	f_data.close();

}

//-----------------------------------------------------------------------------------------------
//                           read in old weights
//-----------------------------------------------------------------------------------------------
void read_weights(NNvariables &NNvar, SIMvariables &SIMvar)
{
	ifstream f_weight;
	int layers, input, output, *arch;
	int activation;
	int check_dim;
	int i,j,k;

	f_weight.open(SIMvar.oldweights.c_str(),ifstream::in);

	if (f_weight.is_open()){
		cout<<"Reading in weights file '"<<SIMvar.oldweights<<"' \n";

		check_dim = 0;
		
		//number of hidden layers
		f_weight>>layers;
		if (layers != NNvar.layers) check_dim = 1;
		arch = new int [layers+2];
		
		// number of inputs 
		f_weight>>input; 			
		if (input != NNvar.input) check_dim = 1;
		arch[0] = input; 
		
		// number of nodes per hidden layers 
		for(i=0; i<layers; i++){
			f_weight>>arch[i+1]; 
			if (arch[i+1]>MAX_NODES-1){
		                cerr << "Fatal Error : cannot use a NN with more than "<<MAX_NODES-1<<" nodes per layers! \n";
				cerr << "Change variable MAX_NODES in header file or change neural network architecture\n\n";
	        	        cerr << "Exiting programme..."<<endl;
				exit(1);
			}
			if (arch[i+1] != NNvar.arch[i+1]) check_dim = 1;
		}


		// number of outputs 
		f_weight>>output;			
		arch[layers+1]=output;
		if (output != NNvar.output) check_dim = 1;

		// activation for output layer
		f_weight>>activation;
		if(activation != NNvar.activation) check_dim = 1;

		// check if NN architecture matches input
		if (check_dim == 1){
			cerr <<"!! ERROR !! Reading weights file: \n";
			cerr <<"NN architecture in input file and in weights file do not match!\n\n";
			cerr <<"input dim: "<<NNvar.input<<" - weights file: "<<input<<endl;
			cerr <<"hidden layers: "<<NNvar.layers<<" - weights file: "<<layers<<endl;
			for(i=1; i<=NNvar.layers; i++){
				cerr <<"Layer "<<i<<": "<<NNvar.arch[i];
				if(i <= layers){
					cerr<<" - weigths file: "<<arch[i];
				}
				cerr<<endl;
			}
			cerr <<"output dim: "<<NNvar.output<<" - weights file: "<<output<<endl;
			cerr <<"activation of output layer: "<<NNvar.activation<<" - weights file: "<<activation<<"  (0=sigmoid, 1=softmax)"<<endl;
			
			cerr<<endl;
			exit(1);
		}

		// input rescaling 
        	for(i=0; i<NNvar.input; i++){
			f_weight>>NNvar.translation[i]; 
			f_weight>>NNvar.stretch[i]; 
		}
	
		// output rescaling 
                for(i=NNvar.input; i<(NNvar.input+NNvar.output); i++){
                        f_weight>>NNvar.translation[i];
                        f_weight>>NNvar.stretch[i];
                }

		// iteration training 
		f_weight>>SIMvar.nprev_iter;

		// read weights 
		for(i=0; i<NNvar.layers+1; i++){
			for(j=0; j<=NNvar.arch[i]; j++){
				for(k=1; k<=NNvar.arch[i+1]; k++){
					NNvar.w[i*MAX_NODES*MAX_NODES+j*MAX_NODES+k]=0;
					f_weight>>NNvar.w[i*MAX_NODES*MAX_NODES+j*MAX_NODES+k];
					//INDICES: this is what Elia uses in nn.cpp because the network is summed up in a different way without using blas!
					//NNvar.w[i*MAX_NODES*MAX_NODES+k*MAX_NODES+j]=0;
					//f_weight>>NNvar.w[i*MAX_NODES*MAX_NODES+k*MAX_NODES+j];
				}
			}
		}


			

	}
	else{
		cerr << "\nFatal Error: cannot open weights file '"<<SIMvar.oldweights<<"' \n";
		cerr << "Exiting programme..."<<endl<<endl;
		exit(1);
	}

	f_weight.close();
}





//-----------------------------------------------------------------------------------------------
//                           reading the input file                               
//-----------------------------------------------------------------------------------------------
void read_input(const string &filename, NNvariables &NNvar, SIMvariables &SIMvar)
{
	ifstream inFile;
	char dummy_char[256];
	inFile.open(filename.c_str(),ifstream::in);

	if (inFile.is_open()){
		//the first 3 lines are comment lines
		inFile.getline(dummy_char,256);
		inFile.getline(dummy_char,256);
		inFile.getline(dummy_char,256);

		//Get rest of the input over a parser
		string string1;
		char instring[201],compstring[201];
		unsigned long linenumber = 0;
		getline(inFile,string1);
		unsigned long run,runi,runstart;
		runstart=16;

		//do this for each line
		do{linenumber++;
			strcpy(instring,string1.c_str());
			char *stringpointer;
			unsigned short readcolumn=0;
			char *col[256];
			strcpy(compstring,instring);
			stringpointer = instring;
			//First 16 characters in col[0]
			//End the string at character no. 16
			instring[15]=0;
			//copy complete string to col[0]
			col[readcolumn]=stringpointer;
			//find the first not empty character to start parsing the rest
			for(runi=16;runi<201;runi++){
				if(instring[runi]!=' '){
					stringpointer = &(instring[runi]);
					runstart=runi;
					break;
				}
				if(runi==200){
					//fpscreen<<"\n\t!ERROR stringpointer could not be set\n\n";
					cerr<<"\n\t!ERROR stringpointer could not be set\n\n";
					stringpointer=&(instring[16]);
					runstart=16;
				}
			}
			//stringpointer = &(instring[16]);  //this does not work if the 17th character is a empty!!
			readcolumn++;
			for(run=runstart;run<201;run++){
				//if you are at the end of the line
				if(instring[run]==0){
					col[readcolumn]=stringpointer;
					stringpointer = &(instring[run+1]);
					readcolumn++;
					break;
				}
				//get values separated by space ' '
				if(instring[run]==' ' && instring[run+1]!=' '){
					instring[run]=0;
					col[readcolumn]=stringpointer;
					stringpointer=&(instring[run+1]);
					readcolumn++;
				}
			}

			// PROCESS THIS LINE (this needs to be done individually!!!)
	
  			// linenumber ... which line is processes
  			// readcolumn ... how many columns does this line have
  			// col[0] bis col[readcolumn-1] ... content of a respective column
  			// if column i contains a number: value is the given by atof(col[i-1]) ... as of type 'double'
			//                                                      atoi(col[i-1]) ... as of type 'int'
			
			if(0==strncmp(compstring,"ran_seed        ",16)){
                        	NNvar.seed=atoi(col[1]);
                        	if(NNvar.seed<=0){
                                	srand(time(NULL));
                                	NNvar.seed=(rand()%1000 +1);
                        	}
                        	continue;
                	}
			if(0==strncmp(compstring,"hidden nodes    ",16)){
				if(NNvar.layers < 0){
					cerr << "Fatal Error : number of 'hidden layers' must be defined in input file before 'hidden nodes'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int nli=0; nli<NNvar.layers;nli++){
						//first one is input layer, so start at nli+1
						NNvar.arch[nli+1] = atoi(col[nli+1]);
						if (NNvar.arch[nli+1]>MAX_NODES-1){
							cerr << "Fatal Error : cannot use a NN with more than "<<MAX_NODES-1<<" nodes per layers! \n";
							cerr << "Change variable MAX_NODES in nn-class-train.h or chnage neural network architecture\n\n";
	        	       				cerr << "Exiting programme..."<<endl;
							exit(1);
						}
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"hidden layers   ",16)){
				NNvar.layers = atoi(col[1]);
				NNvar.arch = new int [NNvar.layers+2];
				continue;
			}
			if(0==strncmp(compstring,"input           ",16)){
				if(NNvar.layers < 0){
					cerr << "Fatal Error : number of 'hidden layers' must be defined in input file before 'input'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					NNvar.arch[0] = atoi(col[1]);
					NNvar.input = atoi(col[1]);
				}
				continue;
			}
			if(0==strncmp(compstring,"output          ",16)){
				if(NNvar.layers < 0){
					cerr << "Fatal Error : number of 'hidden layers' must be defined in input file before 'output'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					NNvar.arch[NNvar.layers+1] = atoi(col[1]);
					NNvar.output = atoi(col[1]);
				}
				continue;
			}
			if(0==strncmp(compstring,"activation      ",16)){
				NNvar.activation = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"run option      ",16)){
				SIMvar.option = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"iterations      ",16)){
				SIMvar.niter = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"save frequency  ",16)){
				SIMvar.nsave = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"optimizer max   ",16)){
				SIMvar.noptmax = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"train input     ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.trainfile = s;
				continue;
			}

			if(0==strncmp(compstring,"test input      ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.testfile = s;
				continue;
			}
			if(0==strncmp(compstring,"test output     ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.testoutf = s;
				continue;
			}
			if(0==strncmp(compstring,"test error      ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.testerrf = s;
				continue;
			}




			if(0==strncmp(compstring,"weights input   ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.oldweights = s;
				continue;
			}
			if(0==strncmp(compstring,"weights output  ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.newweights = s;
				continue;
			}

			if(0==strncmp(compstring,"cost output     ",16)){
				string s = col[1];
				// clean-up white spaces
    				s.erase(remove_if(s.begin(), s.end(), ::isspace),s.end());
				SIMvar.costfile = s;
				continue;
			}

			if(0==strncmp(compstring,"nthreads        ",16)){
				SIMvar.nthreads = atoi(col[1]);
				continue;
			}




		}while (getline(inFile,string1));  //do until end of the file
	}
	else{
		cerr << "Fatal Error: cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl<<endl;
		exit(1);
	}


}


//----------------------------Class NNvariables------------------------------------

//Constructor
NNvariables::NNvariables()
{
	this->seed = 42;
	this->activation = 1; // default is softmax 

	this->layers = -1;
	this->input = -1;
	this->output = -1;
}

//Destructor
NNvariables::~NNvariables(){}


SIMvariables::SIMvariables()
{
	this->option = 0;
	this->niter = 1;
	this->nprev_iter = 0;
	this->nsave = 1;
	this->noptmax = 1;

	this->totalpoints = 0;

	this->oldweights = "weights_old.dat";
	this->newweights = "weights.dat";

	this->testoutf = "output.dat";
	this->testerrf = "error.dat";

	this->costfile = "cost.dat";

	this->nthreads = 1;
}

SIMvariables::~SIMvariables(){}

