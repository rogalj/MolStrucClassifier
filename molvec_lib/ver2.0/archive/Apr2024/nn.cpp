#include "header.h"
#include "dafed.h"


// compute NN for a set of atoms 
/////////////////////////////////////////////////////////////////
void DAFED::get_NN(CAtom *molecules,int center, int &natoms, NNvariables &var){ 
	//cout << "Hey, I'm get_NN" << endl;
	double *inputs   = new double [var.input];
	double *outputs  = new double [var.output];
	double *gradients=new double [var.output*var.input];
	


	for( int kato=0;kato<natoms;kato++){
		if(molecules[kato].eletype !=center)
			continue;
		for( int i=0;i<var.input;i++)inputs[i]=molecules[kato].sfg[i];
		DAFED::NN_network(inputs, outputs, gradients, var);
		for( int i=0;i<var.output;i++)molecules[kato].NNout[i]=outputs[i];	
		for(int i=0;i<(var.input*var.output);i++)molecules[kato].NNgrad[i]=gradients[i];	//BCC gradient [0,12], FCC [13,25], HCP [26,38]
	}

	delete [] inputs; 
	delete [] outputs;
	delete [] gradients; 

}
void DAFED::get_NN_sym_new(CAtom *molecules,int center, int &natoms, NNvariables &var){ 
	//cout << "Hey, I'm get_NN" << endl;
	double *inputs   = new double [var.input];
	double *outputs  = new double [var.output];
	double *gradients=new double [var.output*var.input];

	for( int kato=0;kato<natoms;kato++){
		if(molecules[kato].symtype !=center)
			continue;
		for( int i=0;i<var.input;i++)inputs[i]=molecules[kato].sfg[i];
		DAFED::NN_network(inputs, outputs, gradients, var);
		for( int i=0;i<var.output;i++)molecules[kato].NNout[i]=outputs[i];	
		for(int i=0;i<(var.input*var.output);i++)molecules[kato].NNgrad[i]=gradients[i];	//BCC gradient [0,12], FCC [13,25], HCP [26,38]
	}

	delete [] inputs; 
	delete [] outputs;
	delete [] gradients; 

}
/////////////////////////////////////////////////////////////////
void DAFED::get_NN_norm(CAtom *molecules,int &natoms, NNvariables &var){ 

	double *inputs   = new double [var.input];
	double *outputs  = new double [var.output];
	double *gradients=new double [var.output*var.input];


	for( int kato=0;kato<natoms;kato++){ 

		for( int i=0;i<var.input;i++)inputs[i]=molecules[kato].sfg[i];
		DAFED::NN_network(inputs, outputs, gradients, var);
		for( int i=0;i<var.output;i++)molecules[kato].NNout[i]=outputs[i];	
		for(int i=0;i<(var.input*var.output);i++)molecules[kato].NNgrad[i]=gradients[i];	//BCC gradient [0,12], FCC [13,25], HCP [26,38]

		molecules[kato].NNoutsum = 0.0;
		for(int i=0;i<var.output;i++){
			molecules[kato].NNoutsum += molecules[kato].NNout[i];
		}
		for(int i=0;i<var.output;i++){
			molecules[kato].NNoutnorm[i] = molecules[kato].NNout[i]/molecules[kato].NNoutsum;
		}

	}

	delete [] inputs; 
	delete [] outputs;
	delete [] gradients; 

}




// compute NN for a set of atoms 
/////////////////////////////////////////////////////////////////
void DAFED::get_NN_test(double **sfg,int &natoms, double **strc, double **nngrad, NNvariables &var){ 

	double *inputs   = new double [var.input];
	double *outputs  = new double [var.output];
	double *gradients=new double [var.output*var.input];

	for( int kato=0;kato<natoms;kato++){ 

		for( int i=0;i<var.input;i++)inputs[i]=sfg[kato][i];
		DAFED::NN_network(inputs, outputs, gradients, var);
		for( int i=0;i<var.output;i++)strc[kato][i]=outputs[i];	
		for(int i=0;i<var.input;i++)nngrad[kato][i]=gradients[i];	//BCC gradient [0,12], FCC [13,25], HCP [26,38]
	}

	delete [] inputs; 
	delete [] outputs;
	delete [] gradients; 
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------------------------
// Activation functions 
//-----------------------------------------------------------------------------------------------

// sigmoid
double DAFED::h_class(double x) 
{ 
        double temp;
        temp=(double)1/((double)1+exp(-x));
        return temp;
}

//sigmoid derivative
double DAFED::hprime_class(double x)
{ 
        double temp;
        temp=(double)1/(((double)1+exp(-x))*((double)1+exp(x))); //h_class(x)*((double)1 - h_class(x));
        return temp;
}

// softmax
double DAFED::softmax(double x, double sumexpx)   
{
	double temp;
	temp = exp(x)/sumexpx;
	return temp; 
}

//softmax derivative
double DAFED::softmaxprime(double xi, double xj, int i, int j, double sumexpx)
{
	double tmpi,tmpj;
	double p;
	tmpi = softmax(xi,sumexpx);
	tmpj = softmax(xj,sumexpx);
	if(i==j) p=1.0;
	else p=0.0;
	return (tmpi*(p-tmpj));
}


// Step NN  
///////////////////////////////////////////////////////////////////
void DAFED::NN_step(NNvariables &var,int lay, double *temp){

	for(int i=0; i<MAX_NODES;i++)temp[i]=0;
	
	temp[0]=1.0;

    for(int i=1; i<=var.arch[lay-1] && lay==1; i++)temp[i]=var.nodes[(lay-1)*MAX_NODES+i];
    for(int i=1; i<=var.arch[lay-1] && lay!=1; i++)temp[i]=h_class(var.nodes[(lay-1)*MAX_NODES+i]);

	for(int i=1; i<var.arch[lay]+1; i++){ 
		var.nodes[lay*MAX_NODES+i]=0.0;
		for(int k=0; k<var.arch[lay-1]+1; k++){
			var.nodes[lay*MAX_NODES+i]+=var.w[(lay-1)*MAX_NODES*MAX_NODES+i*MAX_NODES+k]*temp[k];
		}
	}

        var.nodes[lay*MAX_NODES]=1;

}


// Step NN gradient 
/////////////////////////////////////////////////////////////////
void DAFED::NN_stepGradient(NNvariables &var,int lay, double *temp){ 


	//Gradient version
	for(int i=0; i<var.input; i++)
	{
		for(int j=0; j<MAX_NODES; j++)temp[j]=var.Gnodes[(lay-1)*MAX_NODES*MAX_NODES+j*MAX_NODES+i];

        for(int j=1; j<=var.arch[lay-1] && lay!=1; j++)temp[j]=temp[j]*hprime_class(var.nodes[(lay-1)*MAX_NODES+j]);

		for(int j=1; j<var.arch[lay]+1; j++){
			var.Gnodes[lay*MAX_NODES*MAX_NODES+j*MAX_NODES+i]=0.0;	
			for(int k=0; k<var.arch[lay-1]+1; k++){
				var.Gnodes[lay*MAX_NODES*MAX_NODES+j*MAX_NODES+i]+=var.w[(lay-1)*MAX_NODES*MAX_NODES+j*MAX_NODES+k]*temp[k];
			}
		}

        var.Gnodes[lay*MAX_NODES*MAX_NODES+i]=0;
	}
}


// Hidden layers 
/////////////////////////////////////////////////////////////////
void DAFED::NN_hidden(NNvariables &var, int lay){

        double *temp; 
	temp = new double [MAX_NODES];
	
	NN_step(var,lay,temp);
	NN_stepGradient(var,lay,temp);

	delete [] temp;
	if(lay<var.layers+1)DAFED::NN_hidden(var,lay+1);
}


// Forward NN  
/////////////////////////////////////////////////////////////////
void DAFED::NN_network(double *inputs, double *outputs, double *gradient, NNvariables &var){

	double sumexpx = 0;
	double xvalmax = 0; 
	// rescaling and initialize 
    var.nodes[0]=1.0;
    for(int i=1; i<=var.input; i++)var.nodes[i]= (inputs[i-1]-var.translation[i-1])*var.stretch[i-1];
	
	for(int i=0; i<var.input; i++)
	{
		for(int j=0; j<=var.input; j++)var.Gnodes[j*MAX_NODES+i]=0;
		var.Gnodes[(i+1)*MAX_NODES+i]=1.0;
	}

    DAFED::NN_hidden(var,1);

	if (var.activation == 0){
		for(int i=1; i<=var.output; i++)outputs[i-1]=h_class(var.nodes[(var.layers+1)*MAX_NODES+i]);
	}
	else if (var.activation == 1){
		xvalmax = -1e10;
		for(int i=1; i<=var.output; i++) xvalmax = max(xvalmax,var.nodes[(var.layers+1)*MAX_NODES+i]);
		sumexpx = 0.0;
		for(int i=1; i<=var.output; i++) sumexpx += exp(var.nodes[(var.layers+1)*MAX_NODES+i]-xvalmax); 
		for(int i=1; i<=var.output; i++) outputs[i-1]=softmax(var.nodes[(var.layers+1)*MAX_NODES+i]-xvalmax,sumexpx);
	}
	else{
		cerr << "\nError! Something is wrong with variable var.activation in function NN_network (1)...exiting programme\n";
		exit(1);
	}
	
	if(var.activation == 0){
		for(int i=0; i<var.output; i++){
			for(int j=0; j<var.input; j++)gradient[i*var.input+j]=hprime_class(var.nodes[(var.layers+1)*MAX_NODES+i+1])*var.Gnodes[(var.layers+1)*MAX_NODES*MAX_NODES+(i+1)*MAX_NODES+j]*var.stretch[j];
		}
	}
	else if (var.activation == 1){  //xvalmax and sumexpx has been computed above
		for(int i=0; i<var.output; i++){
			for(int j=0; j<var.input; j++){
				gradient[i*var.input+j] = 0.0;
				for(int k=0; k<var.output; k++){
					gradient[i*var.input+j] += softmaxprime(var.nodes[(var.layers+1)*MAX_NODES+i+1]-xvalmax,var.nodes[(var.layers+1)*MAX_NODES+k+1]-xvalmax,i,k,sumexpx)
					*var.Gnodes[(var.layers+1)*MAX_NODES*MAX_NODES+(k+1)*MAX_NODES+j]*var.stretch[j];
				}
			}
		}
	}
	else{
		cerr << "\nError! Something is wrong with variable var.activation in function NN_network (2)...exiting programme\n";
		exit(1);
	}

}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


//read input weight.txy 
/////////////////////////////////////////////////////////////////
void DAFED::read_weight( NNvariables &var)
{
	double dummy;
	ifstream f_weight;

	f_weight.open("weight.txt",ifstream::in); 
	
	if (f_weight.is_open()){
		// number of hidden layers	
		f_weight>>var.layers;			

	
		// Initialiaziation arrays 
		var.arch  = new int [var.layers+2]; 				
		var.w     = new double [(var.layers+1)*MAX_NODES*MAX_NODES];	for(int i=0; i<(var.layers+1)*MAX_NODES*MAX_NODES; i++)var.w[i]=0;
		var.nodes = new double [(var.layers+2)*MAX_NODES];		for(int i=0; i<(var.layers+1)*MAX_NODES; i++)var.nodes[i]=0;
		var.Gnodes= new double [(var.layers+2)*MAX_NODES*MAX_NODES];	for(int i=0; i<(var.layers+1)*MAX_NODES*MAX_NODES; i++)var.Gnodes[i]=0;

		// number of inputs 
		f_weight>>var.input; 			
		var.arch[0] = var.input; 
		
		// number of nodes per hidden layers 
		for(int i=0;i<var.layers;i++){
			f_weight>>var.arch[i+1]; 
			if (var.arch[i+1]>MAX_NODES-1){
		                cerr << "Fatal Error : cannot use a NN with more than "<<MAX_NODES-1<<" nodes per layers! \n";
				cerr << "Change variable MAX_NODES in dafed.h or neural networ\n\n";
	        	        cerr << "Exiting programme..."<<endl;
				exit(1);
			}
		}


		// number of outputs 
		f_weight>>var.output;			
		var.arch[var.layers+1]=var.output;
		// activation for output layer
		f_weight>>var.activation;
		
		var.translation = new double [var.input];
		var.stretch    = new double [var.input];

		// input rescaling 
        	for(int i=0; i<var.input; i++){
			f_weight>>var.translation[i]; 
			f_weight>>var.stretch[i]; 
		}
	
		// output rescaling (no for classification 
                for(int i=0; i<var.output; i++){
                        f_weight>>dummy;
                        f_weight>>dummy;
                }

		// iteration training 
		f_weight>>dummy;
	
		// read weights 
		for(int i=0; i<var.layers+1; i++)
		{
			for(int j=0; j<=var.arch[i]; j++)
			{
				for(int k=1; k<=var.arch[i+1]; k++)
				{
					var.w[i*MAX_NODES*MAX_NODES+k*MAX_NODES+j]=0;
					f_weight>>var.w[i*MAX_NODES*MAX_NODES+k*MAX_NODES+j];
				}
			}
		}


	
	}
	else{
		cerr << "Fatal Error : cannot open the file weight.txt \n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}

	
	f_weight.close();

}


