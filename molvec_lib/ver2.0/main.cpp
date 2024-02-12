///////////////////////////////////////////////////////////////////////////////////////////
//
//
// Determine Structure using Symmetry Functions + Neural Network
// by Daisuke Kuroshima and Jutta Rogal
//
// Daisuke.Kuroshima@nyu.edu
// jutta.rogal@nyu.edu
//
// This code calculates analytical and numerical derivative of the symmetry function.
// The neural network produces a vector that has an entry for each structure.
// 
//
//
///////////////////////////////////////////////////////////////////////////////////////////

#include "header.h"
#include "dafed.h"
#include "readinput.h"
#include <chrono>
using namespace std;


int main(int argc, char *argv[])
{

	int i,j;

	string filename="";
	string filetype="";
	string inputfile="";

	double **x; 		//atom positions
	int *G_ID;
	double box[6]; 		//orthogonal box
	double boxinfo[9]; 		//orthogonal box
	int natoms;		//no of atoms
	int timeframe = 0;	//counter of the step
	int timestep;		//no of timestep
	int *moltype;           //molecule number
	int *eletype;           //element number
	int *symtype;           //element number
	DAFED::CParameter parameter;	//contains some of the parameters for cutoff and symmetry functions
	DAFED::NNvariables NNvar;	//variables for the NN;
	DAFED::EXvariables exvar[2];



	natoms =0;
	timestep = 0;

	for (i=0;i<3;i++){
		box[i] = 0.0;
	}

	if (argc == 4)
	{
		filename = string(argv[1]);
		filetype = string(argv[2]);
		inputfile = string(argv[3]);
	}
	else{
		cout<<"Programme must be called with 3 arguments: filname, filetype, and inputfile!"<<endl;
		cout<<"Exiting programme..."<<endl;
		exit(1);
	}
	if (filetype.compare("lammps") != 0 && filetype.compare("confile") != 0 && filetype.compare("poscar") != 0  && filetype.compare("lammps_aq6") != 0 && filetype.compare("lammps_vec")&& filetype.compare("lammps_vec_tri") && filetype.compare("lammps_sym_new") != 0 && filetype.compare("lammps_sym_new_tri") != 0){
		cout<<"Filetype must be either 'lammps' or 'confile' or 'poscar' or 'lammps_aq6' or 'lammps_vec' or 'lammps_vec_tri' or 'lammps_sym_new' or 'lammps_sym_new_tri'!"<<endl;
		cout<<"Exiting programme..."<<endl;
		exit(1);
	}

	
	get_natoms_lammps(filename,natoms);
	x = new double* [natoms];
	G_ID = new int [natoms];
	moltype = new int [natoms];   //moltype is id of the molecule
        eletype = new int [natoms];   //eletype is element id. 1:H, 2:O, 3:C, 4:N1, 5:N2
        symtype = new int [natoms];   //symmetry function type 0:H, 2:O, 1:C, 3:N1, 4:N2
	for (i=0;i<natoms;i++){
		x[i] = new double [3];
		moltype[i] = 0;
		eletype[i] = 0;
		symtype[i] = 0;
	}
	//check Jutta Rogal's symmetryfunc_basic_deriv directory if you would like to use in atomic system. 
	if (filetype.compare("lammps") == 0 || filetype.compare("lammps_aq6") == 0)
	{
		cout << "This type is no longer accepted. Check the older version";
		cout<<"Exiting programme..."<<endl;
		exit(1);
	}
	//when the system is in orthogonal box
	else if(filetype.compare("lammps_vec") == 0){
		readParticleFile_sa_vec_ortho(filename,x,natoms,moltype,eletype, box, timeframe, timestep, boxinfo, G_ID);
                parameter.triclinic = false;
	}
	//when the system is in triclinic box
	else if(filetype.compare("lammps_vec_tri") == 0){
		readParticleFile_sa_vec_tri(filename,x,natoms,moltype,eletype, box, timeframe,timestep, boxinfo, G_ID);
                parameter.triclinic = true;
	}
	
	//when the system is in orthogonal box
	else if(filetype.compare("lammps_sym_new") == 0){
		readParticleFile_sa_sym_new(filename,x,natoms,moltype,eletype,symtype, box, timeframe, timestep, boxinfo, G_ID);	
                parameter.triclinic = false;
	}
	//when the system is in triclinic box
	else if(filetype.compare("lammps_sym_new_tri") == 0){
		readParticleFile_sa_sym_new_tri(filename,x,natoms,moltype,eletype,symtype, box, timeframe, timestep, boxinfo, G_ID);
		parameter.triclinic = true;
	}

	else
	{
		cout<<"Filetype ERROR'!"<<endl;
		cout<<"Exiting programme..."<<endl;
		exit(1);
	}

	//------------------------------------------------------------------
	//XXXXXX   should be read from a file and then set   XXXXXX
	//parameter.neighbourdistance = 4.0;
	//this is the larger cutoff for the symmetry functions
	////DK: for urea, rmin = 9.8, rmax = 10, for nicotinamide, rmin = 6.8, rmax = 7.0

//	parameter.npairs = 0;

	//a cutoff function with rmin=6.2 and rmax=6.4 gives similar values to Fermi function 

	//for the symmetry function - needs to be read in from file
	//DK:added
	for(i=0;i<parameter.nex;i++){
		exvar[i].dQ = new double* [natoms];
		for(j=0;j<natoms;j++){
			exvar[i].dQ[j] = new double [3];
		}
	}
	if(filetype.compare("lammps_vec") == 0 || filetype.compare("lammps_vec_tri") == 0)
	{
		parameter.sym_new_flag = false;
		read_input(inputfile, parameter);							//read user input file
		DAFED::test_vec(x,natoms,eletype, symtype, moltype, box,parameter,filename,boxinfo, timestep, G_ID);	//inthis function, we run symmetry function calculation
	}

	if(filetype.compare("lammps_sym_new") == 0 || filetype.compare("lammps_sym_new_tri") == 0 )
	{
		read_input(inputfile, parameter);
		parameter.sym_new_flag = true;
		DAFED::test_vec(x,natoms,eletype, symtype, moltype, box,parameter,filename,boxinfo, timestep, G_ID);     //inthis function, we run symmetry function calculation
	}
	return 0;
}

//In this file, we read the input file. Make sure you have all the information
void DAFED::read_input(const string &inputfile, DAFED::CParameter &parameter)
{
        ifstream if_input;

        string string1;
	int current_sf = 0;					//counter of current symmetry function we are calculating 
	int vec_type = 0;					//vec_type starts from 0 up to number of vectors user used 
	int id1, id2;						//first and second element of the vector 	
        if_input.open(inputfile,ifstream::in);                //read INPUT.dat
        if(if_input.is_open()){                                 //checking if the file exist
        	do                                              //check all the line
                {
                        istringstream line(string1);
			if (string1.substr(0,16) =="num_sf          ")	//number of symmetry function
			{
				string test;
				line >> test;
				line >> test;
				parameter.nsfg = stoi(test);
				parameter.symftypeLst = new int[parameter.nsfg];	//vec[0] amd vec[1] of each symf type;
				parameter.pointflag = new bool[parameter.nsfg];
				parameter.g2g3flag = new bool[parameter.nsfg];
	
				parameter.RskappaLst = new double[parameter.nsfg];
				parameter.etaLst = new double[parameter.nsfg];
				for (int i = 0; i < parameter.nsfg;i++)
				{
					parameter.symftypeLst[i] = -1;
					parameter.pointflag[i] = -1;
					parameter.g2g3flag[i] = -1;
					parameter.RskappaLst[i] = -1;
					parameter.etaLst[i] = -1;	
				}
			}
			
			if (string1.substr(0,16) =="num_pvsf_type   ")	//number of type of point vector symmetry function
			{
				string test;
				line >> test;
				line >> test;
				parameter.symftype = stoi(test) + 1;	//need to add 1 for point representation 
				parameter.vectype = new int[2*parameter.symftype];	//vec[0] amd vec[1] of each symf type;
				for (int i = 0; i < 2*parameter.symftype;i++)
					parameter.vectype[i] = -1;
			}	

			if (string1.substr(0,16) =="number_of_para  ")	//type of point vector symmetry function
			{
				string test;
				line >> test;
				line >> test;
				int num = stoi(test);
				getline(if_input,string1);
				istringstream line2(string1);
				if (string1.substr(0,16) =="point           ")	//checking point-vector representation type
				{
					line2 >> test;
					line2 >> test;
					if(test =="true")			//if it is point representation
					{
						for(int i = current_sf;i<current_sf+num; i++)
						{
							parameter.pointflag[i] = true;
							parameter.symftypeLst[i] = 0;
						}
					}
					else					//if it is vector representation
					{
						for(int i = current_sf;i<current_sf+num; i++)
							parameter.pointflag[i] = false;	
						
						getline(if_input,string1);
						istringstream line3(string1);
						
						if (string1.substr(0,16) =="vec_id          ")		//storing vector ID and checking if the set was used before or not
						{
							bool match_flag = false;
							line3 >> test;
							line3 >> test;
							id1 = stoi(test);
							line3 >> test;
							id2 = stoi(test);
							int tmp_vec_type = vec_type;
							for(int i = 2; i < parameter.symftype*2;i+=2)
							{
								if(id1 == parameter.vectype[i] && id2==parameter.vectype[i+1])
								{
									tmp_vec_type = parameter.symftypeLst[2*i];
									match_flag = true;
									break;
								}
							}

							if (tmp_vec_type == vec_type && not(match_flag))
							{
								vec_type++;		//update the number of vec_type
								parameter.vectype[2*vec_type] = id1;
								parameter.vectype[2*vec_type+1] = id2;
								for(int i = current_sf;i<current_sf+num; i++)
									parameter.symftypeLst[i] = vec_type;
							}
							else
								for(int i = current_sf;i<current_sf+num; i++)
									parameter.symftypeLst[i] = tmp_vec_type;
						
						}

						else exit(1);
					}
				}
				else exit(1);
				getline(if_input,string1);
				istringstream line6(string1);
				
				if (string1.substr(0,16) =="G2G3            ")		//checking if it is G2 or G3. Then store parameters
				{
					line6 >> test;
					line6 >> test;
					if(test =="G2")
					{
						for(int i = current_sf;i<current_sf+num; i++)
							parameter.g2g3flag[i] = true;
					
						getline(if_input,string1);
						istringstream line3(string1);
						if (string1.substr(0,16) =="Rs              ")
						{
							line3 >> test;
							for(int i = current_sf;i<current_sf+num; i++)
							{
								line3 >> test;
								parameter.RskappaLst[i] = stod(test);
							}
						}
						else exit(1);
						getline(if_input,string1);
						istringstream line4(string1);

						if (string1.substr(0,16) =="eta             ")
						{
							line4 >> test;
							for(int i = current_sf;i<current_sf+num; i++)
							{
								line4 >> test;
								parameter.etaLst[i] = stod(test);
							}
						}
					}
					else
					{
						for(int i = current_sf;i<current_sf+num; i++)
							parameter.g2g3flag[i] = false;

						getline(if_input,string1);
						istringstream line3(string1);
						if (string1.substr(0,16) =="kappa           ")
						{
							line3 >> test;
							for(int i = current_sf;i<current_sf+num; i++)
                                                        {
                                                                line3 >> test;
                                                                parameter.RskappaLst[i] = stod(test);
								parameter.etaLst[i] = -1;
                                                        }
						}
						else exit(1);
					}
				}
				else exit(1);


				current_sf +=num;
			}

			if (string1.substr(0,16) =="COvectype       ")  //for old code. Checking the first vector
			{
				string test;
				line >> test;                   //separate by space
				
				line >> test;
				parameter.COvectype[0] = stoi(test);
				line >> test;
				parameter.COvectype[1] = stoi(test);
			}
			if (string1.substr(0,16) =="NNvectype       ")  //for old code. Checking the second vector
			{
				string test;
				line >> test;                   //separate by space
				
				line >> test;
				parameter.NNvectype[0] = stoi(test);
				line >> test;
				parameter.NNvectype[1] = stoi(test);	
			}
			if (string1.substr(0,16) =="triclinic       ")  //checking the system box type
			{
				string test;
				line >> test;                   //separate by space
				
				line >> test;
				if ((test=="true" && parameter.triclinic == false) || (test=="false" && parameter.triclinic == true)){
					cout << "File type and INPUT file triclinic parameter is not matching!" << endl;
					exit(1);
				}
			}
			if (string1.substr(0,16) =="center atom     ")  //getting the "center" atom in the molecule
			{
				string test;
				line >> test;                   //separate by space
				line >> test;
				
				line >> test;
				parameter.center = stoi(test);
				parameter.vectype[0] = stoi(test);
				parameter.vectype[1] = stoi(test);
			}
			if (string1.substr(0,16) =="atm in molec    ")  //getting the number of atoms in the molecule
			{
				string test;
				line >> test;                   //separate by space
				line >> test;
				line >> test;
				
				line >> test;
				parameter.natm = stoi(test);
			}
			if (string1.substr(0,16) =="cutoffmin       ")  //getting the low end of the cutoff
			{
				string test;
				line >> test;                   //separate by space
				
				line >> test;
				parameter.rmin0 = stod(test);
			}	
			if (string1.substr(0,16) =="cutoffmax       ")  //getting the high end of the cutoff 
			{
				string test;
				line >> test;                   //separate by space
				
				line >> test;
				parameter.rmax0 = stod(test);
			}	
			if (string1.substr(0,16) =="symf_calc_or_NN ")  //checking if user want to calc symf or NN
                        {
				string test;
                                line >> test;                   //separate by space

                                line >> test;
                                if (test=="symf_calc")
                                	parameter.onlysymcal = 1;
				else if(test=="NN")
					parameter.onlysymcal = 0;
				else
					parameter.onlysymcal = 1;
                        }
			if (string1.substr(0,16) =="NNrdf           ")  //checking if user want to calc symf or NN
			{
				string test;
                                line >> test;                   //separate by space

				line >> test;
				if (test == "true")
				{
					parameter.NNrdf = 1;
					getline(if_input,string1);
                                        istringstream line2(string1);

					if (string1.substr(0,16) =="fullNeigh       ")
					{
						line2 >> test;                   //separate by space
						line2 >> test;
						parameter.surface = stoi(test);
					}
				}
				else
					parameter.NNrdf = 0;
			}
		}while (getline(if_input,string1));  //do until end of the file
	}
	else
		cout <<"ERROR" << endl;
}
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//function definition (usually in corresponding cpp file)


//DK: In here, we calculated symmetry functions only using this directory. (stand alone code)
void DAFED::test_vec(double **pos,int &natoms, int *eletype, int *symtype, int *moltype, double *box, DAFED::CParameter &parameter,string &filename, double * boxinfo, int &timestep, int *G_ID)
{
	CAtom *atoms;			//store all the atom info here
	Cmolpoint *molec;		//store all molecule info here
	DAFED::NNvariables NNvar;	//variables for the NN;
	DAFED::read_weight(NNvar);	//read weight file
	int ti;			
	atoms = new DAFED::CAtom[natoms];	
	molec = new DAFED::Cmolpoint[natoms];
	double **ave_NN = new double* [natoms];
        ofstream symf_file;			
	ofstream NNrdf_file;
        ofstream print_NN;
	ofstream new_traj;
	ofstream run_ave;
	ofstream atminfo;
	symf_file.open("symf_result.txt", ofstream::out);
	
	if (parameter.NNrdf)
	{
		if(!parameter.onlysymcal)print_NN.open("print_NNrdf.dat", ofstream::out);
		NNrdf_file.open("NNrdf_result.dump", ofstream::out);
	}
	if(!parameter.onlysymcal)
	{
        	new_traj.open("new_traj.dump", ofstream::out);
        	run_ave.open("run_ave_traj.dump", ofstream::out);
	}
	
	
	for(ti=0;ti<natoms;ti++){
		atoms[ti].vecpoint = new double* [parameter.symftype];  //each type	
		ave_NN[ti] = new double[parameter.nsfg];
		for (int symf = 0; symf< parameter.symftype;symf++)     //loop for xyz for each symftype
			atoms[ti].vecpoint[symf] = new double[3];
		for (int nsfg = 0; nsfg< parameter.nsfg; nsfg++)
			ave_NN[ti][nsfg] = 0;
        }

	for(int mi = 0; mi<natoms;mi++)
	{
		molec[mi].atom_ID = new int[parameter.natm];
		molec[mi].vector_ID = new int[parameter.symftype*2];	//type * 2 
		molec[mi].vecpoint = new double* [parameter.symftype];	//each type
		for (int symf = 0; symf< parameter.symftype;symf++)	//loop for xyz for each symftype
			molec[mi].vecpoint[symf] = new double[3];	
	}
	for(int i=0;i<natoms;i++){
		atoms[i].sfg = new double[parameter.nsfg];
		atoms[i].NNout = new double[NNvar.output];
		atoms[i].NNgrad = new double[NNvar.output*NNvar.input];
	}
	int timeframe = 1;
	//here we read dump file and calculate symf or NN and output to symf_result.txt and/or new_traj.dump
	do
	{
		
	for(ti=0;ti<natoms;ti++){
                for(int j=0;j<3;j++){
	                atoms[G_ID[ti]].pos[j] = pos[ti][j];
		}
		atoms[G_ID[ti]].eletype = eletype[ti];
		atoms[G_ID[ti]].moltype = moltype[ti];
		if(parameter.sym_new_flag)atoms[G_ID[ti]].symtype = symtype[ti];
	}
        atminfo.open("atominfo.txt", ofstream::out);
        for(ti=0;ti<natoms;ti++){					//DK: we get all the atominfo to atominfo.txt. Checking if all information successfully read or not
		for (int k = 0; k<3; k++)
		{
			atminfo << "atom " << ti;
                	atminfo << " moltype is :" << atoms[G_ID[ti]].moltype;
                	atminfo << " eletype is :" << atoms[G_ID[ti]].eletype;
                	if(parameter.sym_new_flag)atminfo << " symtype is :" << atoms[G_ID[ti]].symtype;
			atminfo << " position  is " << atoms[G_ID[ti]].pos[k]<<endl;
		}
	}
	atminfo.close();
		
	//here we calculate neighbor information
	//if you use the new method, sym_new_flag is true, otherwise go to the original
	if(parameter.sym_new_flag)get_allNeighbourDistances_sa_sym_new(atoms,molec,natoms,box,parameter);
	else get_allNeighbourDistances_gen_vec_sa(atoms,molec,natoms,box,parameter);

	//Here we calculate the symmetry function 
	get_vec_symmetry_NN_vec_sa_short_test(atoms,molec,natoms,parameter,NNvar,symf_file);
	
	//Printing out to the new_traj.dump file with NNout results with dump format. 
	for (int mi = 0; mi < parameter.nmol;mi++)
	{
		for (int ti = 1;ti < parameter.natm;ti++)
		{
			if(molec[mi].atom_ID[ti]==-1)break;                   //If it does not exist, we skip it.
			for (int j =0;j < NNvar.output;j++)
				atoms[molec[mi].atom_ID[ti]].NNout[j] = atoms[molec[mi].atom_ID[0]].NNout[j];
		}
	}
	//here we print out to new_traj.dump file
	
	print_traj(atoms, molec, box, natoms, ave_NN, parameter, NNvar, boxinfo, timestep, timeframe, G_ID, new_traj, run_ave, print_NN, NNrdf_file);
	
	if(parameter.sym_new_flag)
	{
		if(parameter.triclinic) readParticleFile_sa_sym_new_tri(filename,pos,natoms,moltype,eletype,symtype,box, timeframe, timestep, boxinfo, G_ID);
		else readParticleFile_sa_sym_new(filename,pos,natoms,moltype,eletype,symtype,box, timeframe, timestep, boxinfo, G_ID);
	}
	else
	{
		if(parameter.triclinic) readParticleFile_sa_vec_tri(filename,pos,natoms,moltype,eletype,box, timeframe, timestep, boxinfo, G_ID);
		else readParticleFile_sa_vec_ortho(filename,pos,natoms,moltype,eletype,box, timeframe, timestep, boxinfo, G_ID);
	}
	}while(timeframe!=-1);

	if (parameter.NNrdf)
	{
		if(!parameter.onlysymcal)print_NN.close();
		NNrdf_file.close();
	}
	if(!parameter.onlysymcal)
	{
        	new_traj.close();
        	run_ave.close();
	}
	symf_file.close();
	exit(0);
}

void DAFED::print_traj(CAtom *atoms, Cmolpoint *molec, double *box, int &natoms, double **ave_NN, CParameter &parameter, NNvariables &NNvar, double *boxinfo, int &timestep, int &timeframe, int *G_ID, ofstream &new_traj, ofstream &run_ave, ofstream &print_NN, ofstream &NNrdf_file)
{	
	if (parameter.NNrdf)
	{	
		if(!parameter.onlysymcal)
		{
			new_traj << "ITEM: TIMESTEP\n" << timestep << endl;
			run_ave << "ITEM: TIMESTEP\n" << timestep << endl;
			new_traj << "ITEM: NUMBER OF ATOMS\n" << natoms << endl;
			run_ave << "ITEM: NUMBER OF ATOMS\n" << natoms << endl;
        		NNrdf_file << "ITEM: TIMESTEP\n" << timestep << endl;
			NNrdf_file << "ITEM: NUMBER OF ATOMS\n" << natoms << endl;

			if(parameter.triclinic)
			{
				new_traj<<"ITEM: BOX BOUNDS xy xz yz pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<" "<<boxinfo[2]<< endl 
					<<boxinfo[3]<<" "<<boxinfo[4]<<" "<<boxinfo[5]<<endl<<boxinfo[6]<<" "<<boxinfo[7]<<" "<<boxinfo[8]<<endl;
				run_ave<<"ITEM: BOX BOUNDS xy xz yz pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<" "<<boxinfo[2]<< endl 
					<<boxinfo[3]<<" "<<boxinfo[4]<<" "<<boxinfo[5]<<endl<<boxinfo[6]<<" "<<boxinfo[7]<<" "<<boxinfo[8]<<endl;
				NNrdf_file << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<" "<<boxinfo[2]<< endl 
					<<boxinfo[3]<<" "<<boxinfo[4]<<" "<<boxinfo[5]<<endl<<boxinfo[6]<<" "<<boxinfo[7]<<" "<<boxinfo[8]<<endl;
			}
			else
			{
				new_traj<<"ITEM: BOX BOUNDS pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<endl<<boxinfo[2]<<" "<<boxinfo[3]<<endl<<boxinfo[4]<<" "<<boxinfo[5]<<endl;
				run_ave<<"ITEM: BOX BOUNDS pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<endl<<boxinfo[2]<<" "<<boxinfo[3]<<endl<<boxinfo[4]<<" "<<boxinfo[5]<<endl;
				NNrdf_file << "ITEM: BOX BOUNDS pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<endl<<boxinfo[2]<<" "<<boxinfo[3]<<endl<<boxinfo[4]<<" "<<boxinfo[5]<<endl;
			}
			new_traj << "ITEM: ATOMS id mol element x y z n_nei ";
			run_ave << "ITEM: ATOMS id mol element x y z n_nei ";
			NNrdf_file <<"ITEM: ATOMS id mol element x y z n_nei dist sec_mol ";
			
			for (int j =0;j < NNvar.output;j++)
			{
				new_traj << "NNout" + to_string(j+1) + " ";
				run_ave << "NNout" + to_string(j+1) + " ";
				NNrdf_file << "NNout" + to_string(j+1) + " ";
			}
			
			new_traj << endl;
			run_ave << endl;
			NNrdf_file << endl;
			int g_id = -1;
	
			for (int i = 0; i < natoms; i++)
			{
				g_id = G_ID[i];
				new_traj << g_id+1 << " " << atoms[g_id].moltype +1 << " " << atoms[g_id].eletype  << " " 
					<< atoms[g_id].pos[0] << " " << atoms[g_id].pos[1] << " " << atoms[g_id].pos[2] << " " << atoms[g_id].n_neighbors << " ";
				run_ave << g_id+1 << " " << atoms[g_id].moltype +1 << " " << atoms[g_id].eletype  << " " 
					<< atoms[g_id].pos[0] << " " << atoms[g_id].pos[1] << " " << atoms[g_id].pos[2] << " " << atoms[g_id].n_neighbors << " ";
				get_NNrdf(atoms, molec, natoms, parameter, box);
				NNrdf_file << g_id+1 << " " << atoms[g_id].moltype +1 << " " << atoms[g_id].eletype  
					<< " " << atoms[g_id].pos[0] << " " << atoms[g_id].pos[1] << " " << atoms[g_id].pos[2] 
					<< " " << atoms[g_id].n_neighbors << " "<< atoms[g_id].d_surface << " " << atoms[g_id].second_mol << " ";

				for (int j =0;j < NNvar.output;j++)
				{
					new_traj << atoms[g_id].NNout[j] << " ";
					ave_NN[g_id][j] += atoms[g_id].NNout[j];
					run_ave << ave_NN[g_id][j]/timeframe <<" ";
					NNrdf_file << ave_NN[g_id][j]/timeframe <<" ";
				}
				new_traj << endl;
				run_ave << endl;
				NNrdf_file << endl;
			}
		//	print_NN << "# " << timestep << endl;
		//	get_print_NNrdf(atoms, molec, natoms, parameter, NNvar, print_NN);				
		}
		else
		{
        		NNrdf_file << "ITEM: TIMESTEP\n" << timestep << endl;
			NNrdf_file << "ITEM: NUMBER OF ATOMS\n" << natoms << endl;
			if(parameter.triclinic)
			{
				NNrdf_file << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<" "<<boxinfo[2]<< endl <<boxinfo[3]
					<<" "<<boxinfo[4]<<" "<<boxinfo[5]<<endl<<boxinfo[6]<<" "<<boxinfo[7]<<" "<<boxinfo[8]<<endl;
			}
			else
			{
				NNrdf_file << "ITEM: BOX BOUNDS pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<endl<<boxinfo[2]<<" "<<boxinfo[3]<<endl<<boxinfo[4]<<" "<<boxinfo[5]<<endl;
			}
			NNrdf_file <<"ITEM: ATOMS id mol element x y z n_nei dist sec_mol ";
			NNrdf_file << endl;
			int g_id = -1;
			for (int i = 0; i < natoms; i++)
			{
				g_id = G_ID[i];
				get_NNrdf(atoms, molec, natoms, parameter, box);
				NNrdf_file << g_id+1 << " " << atoms[g_id].moltype +1 << " " << atoms[g_id].eletype  
					<< " " << atoms[g_id].pos[0] << " " << atoms[g_id].pos[1] << " " << atoms[g_id].pos[2] << " " 
					<< atoms[g_id].n_neighbors << " "<< atoms[g_id].d_surface << " " << atoms[g_id].second_mol << " ";
				NNrdf_file << endl;
			}
		}
	}
	else
	{
		if(!parameter.onlysymcal)
                {
			new_traj << "ITEM: TIMESTEP\n" << timestep << endl;
			run_ave << "ITEM: TIMESTEP\n" << timestep << endl;
			new_traj << "ITEM: NUMBER OF ATOMS\n" << natoms << endl;
			run_ave << "ITEM: NUMBER OF ATOMS\n" << natoms << endl;

			if(parameter.triclinic)
			{
				new_traj<<"ITEM: BOX BOUNDS xy xz yz pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<" "<<boxinfo[2]<< endl 
					<<boxinfo[3]<<" "<<boxinfo[4]<<" "<<boxinfo[5]<<endl<<boxinfo[6]<<" "<<boxinfo[7]<<" "<<boxinfo[8]<<endl;
				run_ave<<"ITEM: BOX BOUNDS xy xz yz pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<" "<<boxinfo[2]<< endl 
					<<boxinfo[3]<<" "<<boxinfo[4]<<" "<<boxinfo[5]<<endl<<boxinfo[6]<<" "<<boxinfo[7]<<" "<<boxinfo[8]<<endl;
			}
			else
			{
				new_traj<<"ITEM: BOX BOUNDS pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<endl<<boxinfo[2]<<" "<<boxinfo[3]<<endl<<boxinfo[4]<<" "<<boxinfo[5]<<endl;
				run_ave<<"ITEM: BOX BOUNDS pp pp pp\n"<<boxinfo[0]<<" "<<boxinfo[1]<<endl<<boxinfo[2]<<" "<<boxinfo[3]<<endl<<boxinfo[4]<<" "<<boxinfo[5]<<endl;
			}
			new_traj << "ITEM: ATOMS id mol element x y z n_nei ";
			run_ave << "ITEM: ATOMS id mol element x y z n_nei ";
			
			for (int j =0;j < NNvar.output;j++)
			{
				new_traj << "NNout" + to_string(j+1) + " ";
				run_ave << "NNout" + to_string(j+1) + " ";
			}
			
			new_traj << endl;
			run_ave << endl;
			int g_id = -1;
	
			for (int i = 0; i < natoms; i++)
			{
				g_id = G_ID[i];
				new_traj << g_id+1 << " " << atoms[g_id].moltype +1 << " " << atoms[g_id].eletype  << " " 
					<< atoms[g_id].pos[0] << " " << atoms[g_id].pos[1] << " " << atoms[g_id].pos[2] << " " << atoms[g_id].n_neighbors << " ";
				run_ave << g_id+1 << " " << atoms[g_id].moltype +1 << " " << atoms[g_id].eletype  << " " 
					<< atoms[g_id].pos[0] << " " << atoms[g_id].pos[1] << " " << atoms[g_id].pos[2] << " " << atoms[g_id].n_neighbors << " ";

				for (int j =0;j < NNvar.output;j++)
				{
					new_traj << atoms[g_id].NNout[j] << " ";
					ave_NN[g_id][j] += atoms[g_id].NNout[j];
					run_ave << ave_NN[g_id][j]/timeframe <<" ";
				}
				new_traj << endl;
				run_ave << endl;
			}
		}
	}
}

//----------------------------------------------------------------------------------------
// read DAFED input file (only for lammps)
void DAFED::read_dafed_input(const string &filename, CParameter &para, EXvariables *exvar,  string &cfile, int &cv_stride, string &dcvfile, int &dcv_stride)
{
	ifstream inFile;
	char dummy_char[256];
	int subtotal = 0;

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

			if(0==strncmp(compstring,"ex_kappa        ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_kappa'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].kappa = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_tau          ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_tau'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].tau = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_gamma        ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_gamma'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].gamma = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_temp         ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_temp'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].temp = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_mass         ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_mass'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].m = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_nrespa       ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_nrespa'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].ggmt.n_respa_ggmt = atoi(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_histo_nbin   ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_histo_nbin'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].histo.nbin = atoi(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_histo_min    ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_histo_min'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].histo.min = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_histo_max    ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_histo_max'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].histo.max = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_lwall        ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_lwall'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].lwall = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_lwall_n      ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_lwall_n'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].lwall_n = atoi(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_lwall_eps    ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_lwall_eps'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].lwall_eps = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_lwall_k      ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_lwall_k'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].lwall_k = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_lwall_cv     ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_lwall_cv'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].lwall_cv = atoi(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_uwall        ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_uwall'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].uwall = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_uwall_n      ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_uwall_n'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].uwall_n = atoi(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_uwall_eps    ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_uwall_eps'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].uwall_eps = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_uwall_k      ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_uwall_k'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].uwall_k = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_uwall_cv     ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_uwall_cv'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].uwall_cv = atoi(col[exi+1]);
					}
				}
				continue;
			}
			//JR start: included for ramping s up/down
			if(0==strncmp(compstring,"ex_ramp_smin    ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_ramp_smin'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].ramp.smin = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_ramp_smax    ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_ramp_smax'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].ramp.smax = atof(col[exi+1]);
					}
				}
				continue;
			}
			if(0==strncmp(compstring,"ex_ramp_ds      ",16)){
				if(para.nex < 0){
					cerr << "Fatal Error : number of extended variables 'ex_number' must be defined in input file before 'ex_ramp_ds'"<<endl;
					cerr << "Change order in input file "<< filename<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				else{
					for(int exi=0; exi<para.nex;exi++){
						exvar[exi].ramp.ds = atof(col[exi+1]);
					}
				}
				continue;
			}
			//JR end: included for ramping s up/down

			if(0==strncmp(compstring,"ex_number       ",16)){
				para.nex = atoi(col[1]);
				if(para.nex > 2){
					cerr << "Fatal Error : too many extended variables, maximum number is 2!"<<endl;
					cerr << "Exiting programme..."<<endl;
					exit(1);
				}
				continue;
			}

			if(0==strncmp(compstring,"ex_seed         ",16)){
                        	para.ex_seed=atoi(col[1]);
                        	//fpscreen<<"Seed read from input \t = "<<pseed<<endl;
                        	if(para.ex_seed>=0){
                                	srand(time(NULL));
                                	para.ex_seed=-(rand()%1000 +1);
                        	}
                        	//fpscreen<<"Seed used \t\t = "<<pseed<<"\n";
                        	continue;
                	}
			if(0==strncmp(compstring,"use_NN          ",16)){
				para.useNN = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"choose_exvar    ",16)){
				para.choose_exvar = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"use_metadyn     ",16)){
				para.useMetadyn = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"gauss_freq      ",16)){
				para.gauss_freq = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"gauss_height    ",16)){
				para.gauss_h = atof(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"gauss_width     ",16)){
				para.gauss_sigma = atof(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"bias_stride     ",16)){
				para.bias_stride = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"read_bias       ",16)){
				para.bias_read = atoi(col[1]);
				para.bias_readfile = col[2];
				para.bias_readfile.erase(remove(para.bias_readfile.begin(),para.bias_readfile.end(),' '),para.bias_readfile.end());
				continue;
			}
			if(0==strncmp(compstring,"CV_file         ",16)){
				cfile=col[1];
				cfile.erase(remove(cfile.begin(),cfile.end(),' '),cfile.end());
				continue;
			}
			if(0==strncmp(compstring,"CV_stride       ",16)){
				cv_stride = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"dCVdr_file      ",16)){
				dcvfile=col[1];
				dcvfile.erase(remove(dcvfile.begin(),dcvfile.end(),' '),dcvfile.end());
				continue;
			}
			if(0==strncmp(compstring,"dCVdr_stride    ",16)){
				dcv_stride = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"lower_wall      ",16)){
				para.lwall = atof(col[1]);
				para.lwall_cv = atoi(col[2]);
				continue;
			}
			if(0==strncmp(compstring,"lwall_para      ",16)){
				para.lwall_n = atoi(col[1]);
				para.lwall_eps = atof(col[2]);
				para.lwall_k = atof(col[3]);
				continue;
			}
			if(0==strncmp(compstring,"use_restraint   ",16)){
				para.use_restraint = atoi(col[1]);
				para.restraint_value = atof(col[2]);
				continue;
			}
			if(0==strncmp(compstring,"restraint_para  ",16)){
				para.restraint_n = atoi(col[1]);
				para.restraint_eps = atof(col[2]);
				para.restraint_k = atof(col[3]);
				continue;
			}

			//DKadded:
			if(0==strncmp(compstring,"G2_point        ",16)){
				para.nsfg2point = atoi(col[1]);
				cout << "para.nsfg2point = " << para.nsfg2point << endl;
				continue;
			}
			if(0==strncmp(compstring,"G2p_Rs          ",16)){
				for (int sf = subtotal; sf < subtotal + para.nsfg2point; sf++)
					para.RskappaLst[sf] = atof(col[sf-subtotal+1]);
				continue;
			}
			if(0==strncmp(compstring,"G2p_eta         ",16)){
				for (int sf = subtotal; sf < subtotal + para.nsfg2point; sf++)
					para.etaLst[sf] = atof(col[sf-subtotal+1]);
				continue;
			}

			if(0==strncmp(compstring,"G2_CO           ",16)){
				para.nsfg2CO = atoi(col[1]);
				subtotal += para.nsfg2point;
				continue;
			}
			if(0==strncmp(compstring,"G2CO_Rs         ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg2CO; sf++)
					para.RskappaLst[sf] = atof(col[sf-subtotal+1]);
				continue;
			}
			if(0==strncmp(compstring,"G2CO_eta        ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg2CO; sf++)
					para.etaLst[sf] = atof(col[sf-subtotal+1]);
				continue;
			}

			if(0==strncmp(compstring,"G2_NN           ",16)){
				para.nsfg2NN = atoi(col[1]);
				subtotal += para.nsfg2CO;
				continue;
			}
			if(0==strncmp(compstring,"G2NN_Rs         ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg2NN; sf++)
					para.RskappaLst[sf] = atof(col[sf-subtotal+1]);
				continue;
			}
			if(0==strncmp(compstring,"G2NN_eta        ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg2NN; sf++)
					para.etaLst[sf] = atof(col[sf-subtotal+1]);
				continue;
			}

			if(0==strncmp(compstring,"G3_point        ",16)){
				para.nsfg3point = atoi(col[1]);
				subtotal += para.nsfg2NN;
				continue;
			}
			if(0==strncmp(compstring,"G3p_Rs          ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg3point; sf++)
				{
					para.RskappaLst[sf] = atof(col[sf-subtotal+1]);
					para.etaLst[sf] = -1;
				}
				continue;
			}

			if(0==strncmp(compstring,"G3_CO           ",16)){
				para.nsfg3CO = atoi(col[1]);
				subtotal += para.nsfg3point;
				continue;
			}
			if(0==strncmp(compstring,"G3CO_Rs         ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg3CO; sf++)
				{
					para.RskappaLst[sf] = atof(col[sf-subtotal+1]);
					para.etaLst[sf] = -1;
				}
				continue;
			}	
			
			if(0==strncmp(compstring,"G3_NN           ",16)){
				para.nsfg3NN = atoi(col[1]);
				subtotal += para.nsfg3CO;
				continue;
			}
			if(0==strncmp(compstring,"G3NN_Rs         ",16)){
				for (int sf = subtotal; sf < subtotal+para.nsfg3NN; sf++)
				{
					para.RskappaLst[sf] = atof(col[sf-subtotal+1]);
					para.etaLst[sf] = -1;
				}
				continue;
			}
			
			
			if(0==strncmp(compstring,"center_atom     ",16)){
				para.center = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"atm_in_molec    ",16)){
				para.natm = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"num_of_NN_out    ",16)){
				para.nnout = atoi(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"eleOfvec1        ",16)){
				para.COvectype[0] = atoi(col[1]);
				para.COvectype[1] = atoi(col[2]);
				continue;
			}
			if(0==strncmp(compstring,"eleOfvec2        ",16)){
				para.NNvectype[0] = atoi(col[1]);
				para.NNvectype[1] = atoi(col[2]);
				continue;
			}
			if(0==strncmp(compstring,"cutoff           ",16)){
				para.rmax0 = atof(col[1]);
				continue;
			}
			if(0==strncmp(compstring,"cutoff_diff      ",16)){
				para.rmin0 = para.rmax0 - atof(col[1]);
				continue;
			}

			//print out text string
			if(0==strncmp(compstring,"#",1))
			{
				//fpscreen<<"\n"<<compstring<<"\n";
			}
		}while (getline(inFile,string1));  //do until end of the file

	}
	else{
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}



//----------------------------------------------------------------------------------------
// random number generator, uniform random number
double DAFED::ran3(int idum)
{ 
	static int inext, inextp;
  	static int iff=0;
  	const int MBIG=1000000000,MSEED=161803398,MZ=0;
  	const double FAC=(1.0/MBIG);
	//  static valarray<int> ma(56);
  	static int ma[56];
  	int i, ii, k, mj, mk;
  	if(idum<0 || iff==0){ 
		iff=1;
    		mj=labs(MSEED-labs(idum));
    		mj%=MBIG;
    		ma[55]=mj;
    		mk=1;
    		for (i=1;i<55;i++){ 
			ii = (21*i) % 55;
      			ma[ii] = mk;
      			mk = mj - mk;
      			if (mk<int(MZ)) mk+=MBIG;
      			mj = ma[ii];
		}
    		for(k=1;k<=4;k++){ 
			for(i=1;i<=55;i++){ 
				ma[i] = ma[i] - ma[1 + ((i+30)%55)];
        			if (ma[i]<int(MZ)) ma[i]+=MBIG;
			}
		}
    		inext=0;
    		inextp=31;
    		//idum=1;		//this must be in here, if you don't have ran3(int idum=1) in  header file
    					//and then the reference must be passed 
    					//and function must always be called with ran3(idum)
	}
  	if (++inext == 56)   inext  =1;
  	if (++inextp ==56)  inextp =1;
  	mj = ma[inext] - ma[inextp];
  	if (mj<int(MZ)) mj+=MBIG;
  	ma[inext]=mj;
  	return  mj *FAC;
}

// Gaussian random number --------------
double DAFED::gauss()
{
        double v1,v2,rsq,fac;
	double grand;
	static int getnew = 0;
	static double savegrand = 0.0;


	if(getnew == 0){
        	do
        	{
                	v1=2.0*ran3()-1.0;
                	v2=2.0*ran3()-1.0;
                	rsq= v1*v1 + v2*v2;
        	}while((rsq>=1) || (rsq==0));
        	fac = sqrt(-2.0*log(rsq)/rsq);

        	//v1*fac and v2*fac are gaussian distributed random numbers (cf. numerical recipies etc.)
		grand = v1*fac;
		savegrand = v2*fac;
		getnew = 1;
	}
	else{
		grand = savegrand;
		getnew = 0;
	}



        return grand;
}


void DAFED::get_Q_vec_lmpneigh_vec_new(CAtom *molecules, Cmolpoint *molec, CNeighvariables &lmpneigh, EXvariables *exvar, CParameter &parameter, double *Qlocal)
{
	int i;
	int ti;

	int natoms = lmpneigh.inum;	//no of owned atoms on this proc (sort of, but is the no for which a neighbour list exists)
	
	int nmol=parameter.nmol;
	int iex;
	int nex=parameter.nex;
	int outNN_index[2];
	int outi;
	
	exvar[0].Q = 0.0;
	exvar[1].Q = 0.0;

	//this function determines sfg and NN and NNgrad
	//for parallel in lammps this has to be called first, then the info on the ghost atoms needs
	//to be exchanged, then we can continue with this function
	//get_symmetry_NN_3_lmpneighbour(molecules,lmpneigh,box,parameter,sf2Rs,sf2eta,sf3kappa,NNvar);

	//this is for nex=2 with 0=BCC and 1=A15	
	//and     for nex=1 with 0=A15
	//
	get_total_derivatives_molvec_vec_lmp(molecules, molec, natoms, parameter, exvar,lmpneigh);
	//determine values of extended variables
	//if this is done parallel in lammps, then this is only the local sum
	for(iex=0;iex<nex;iex++){
		exvar[iex].Qlocal = 0.0;
	}
	for(i=0;i<parameter.nnout;i++){
		Qlocal[i] = 0.0;
	}
	

	if(nex==2){				//this is for nex=2 and  0=BCC and 1=A15
		outNN_index[0] = 0;		//U1
		outNN_index[1] = 1;		//U4
	}
	else if(nex==1){			//this is for nex=1 and 0=U1
		outNN_index[0] = 0;		//A15
	}
	else{
		cerr<<"\n\n\t!!ERROR!! In get_Q_3_lmpneighbour, no. of extended variable nex must be either 1 or 2!\n";
		cerr<<"\tnex = "<<nex<<endl;
		cerr<<"\tExiting programme!\n\n";
		exit(1);
	}

	for(ti=0;ti<natoms;ti++){
		if(molecules[ti].eletype!=parameter.center)continue;
		for(iex=0;iex<nex;iex++){
			outi = outNN_index[iex];
			exvar[iex].Qlocal += molecules[ti].NNout[outi];
		}
		for(i=0;i<parameter.nnout;i++){
			Qlocal[i] += molecules[ti].NNout[i];
		}
	}

	for(iex=0;iex<nex;iex++){
		exvar[iex].Qlocal = exvar[iex].Qlocal/nmol;		//devided by total number of mol in the simulation
	}
	for(i=0;i<parameter.nnout;i++){
		Qlocal[i] = Qlocal[i]/nmol;
	}
	//print value of order parameter
}







//This is one of the old code that accept only two vectors. Calculates symmetry function then run NN. lammps only
void DAFED::get_vec_symmetry_NN_3_lmpneighbour_new(CAtom *atoms, Cmolpoint *molec, CNeighvariables &lmpneigh, double *box, CParameter &parameter, double *Rskappa, double *eta,  NNvariables &NNvar)
{
	int i;
	int nsfg,nsfg2CO,nsfg2NN,nsfg3CO,nsfg3NN,nsfg2point,nsfg3point;
	int natoms = lmpneigh.inum;
        int COvectype[2];
        int NNvectype[2];
	COvectype[0] = parameter.COvectype[0];
	COvectype[1] = parameter.COvectype[1];
	NNvectype[0] = parameter.NNvectype[0];
	NNvectype[1] = parameter.NNvectype[1];

	nsfg=parameter.nsfg;
	nsfg2point=parameter.nsfg2point;
	nsfg2CO=parameter.nsfg2CO;
	nsfg2NN=parameter.nsfg2NN;
        nsfg3point=parameter.nsfg3point;
	nsfg3CO=parameter.nsfg3CO;
	nsfg3NN=parameter.nsfg3NN;

	int nmol = lmpneigh.lmpmol;
	
	for (i = 0; i < nsfg; i++){
		parameter.Rskappa=Rskappa[i];
		parameter.eta=eta[i];
		if (i < nsfg2point){
			symfg2_lmp(atoms, molec, natoms, parameter,lmpneigh);
			for(int mi = 0; mi < nmol;mi++)
				atoms[molec[mi].center].sfg[i] = molec[mi].point2;	
		}
		else if (i < nsfg2point+nsfg2CO){
			symmetryfunc_g2_fc0_vecCO_lmp(atoms, molec, natoms, parameter, COvectype,lmpneigh);
			for(int mi = 0; mi < nmol;mi++)
				atoms[molec[mi].center].sfg[i] = molec[mi].G2vecCO;
		}

		else if (i < nsfg2point+nsfg2CO+nsfg2NN){
			symmetryfunc_g2_fc0_vecNN_lmp(atoms, molec, natoms, parameter, NNvectype,lmpneigh);
			for(int mi = 0; mi < nmol;mi++)
				atoms[molec[mi].center].sfg[i] = molec[mi].G2vecNN;
		}
	
		else if (i < nsfg2point+nsfg2CO+nsfg2NN+nsfg3point){
			symfg3_lmp(atoms, molec, natoms, parameter,lmpneigh);
			for(int mi = 0; mi < nmol;mi++)
				atoms[molec[mi].center].sfg[i] = molec[mi].point3;
		}
		
		else if (i < nsfg2point+nsfg2CO+nsfg2NN+nsfg3point+nsfg3CO){
			symmetryfunc_g3_fc0_vecCO_lmp(atoms, molec, natoms, parameter, COvectype,lmpneigh);
			for(int mi = 0; mi < nmol;mi++)
				atoms[molec[mi].center].sfg[i] = molec[mi].G3vecCO;
		}	
		
		else if (i < nsfg2point+nsfg2CO+nsfg2NN+nsfg3point+nsfg3CO+nsfg3NN){
			symmetryfunc_g3_fc0_vecNN_lmp(atoms, molec, natoms, parameter, NNvectype,lmpneigh);
			for(int mi = 0; mi < nmol;mi++)
				atoms[molec[mi].center].sfg[i] = molec[mi].G3vecNN;
		}
		
		else
			cout << "Check your nsfg and length of Rskappa list!" << endl;
	}
	//calculate NN output and derivatives
	get_NN(atoms,parameter.center,natoms,NNvar);
}

//Calculates symmetry function then run NN. lammps only
void DAFED::get_vec_symmetry_NN_3_lmpneighbour_vec_lmp(CAtom *atoms, Cmolpoint *molec, CNeighvariables &lmpneigh, CParameter &parameter, NNvariables &NNvar)
{
        int i;
        int nsfg;
        nsfg=parameter.nsfg;
	int natoms = lmpneigh.inum;
	int nmol = lmpneigh.lmpmol;
	for (i = 0; i < nsfg; i++){
                symf_multi_vec_lmp(atoms, molec, parameter, lmpneigh, parameter.symftypeLst[i], i);
		for(int mi = 0; mi < nmol;mi++)
                        atoms[molec[mi].center].sfg[i] = molec[mi].symf;
        }
	//calculate NN output and derivatives
	get_NN(atoms,parameter.center,natoms,NNvar);
}

//Calculates symmetry function and NN lammps only
void DAFED::get_vec_symmetry_NN_vec_lmp_short_test(CAtom *atoms, Cmolpoint *molec, CNeighvariables &lmpneigh, CParameter &parameter, NNvariables &NNvar)
{
	int natoms = lmpneigh.inum;
	symf_multi_vec_lmp_short_test(atoms,molec,parameter,lmpneigh);
	
	//calculate NN output and derivatives
	get_NN(atoms,parameter.center,natoms,NNvar);
}



//calculate symf then NN. stand alone only
void DAFED::get_vec_symmetry_NN_vec_sa(CAtom *atoms, Cmolpoint *molec,int &natoms, CParameter &parameter, NNvariables &NNvar)
{
	int i;
	int nsfg;
	nsfg=parameter.nsfg;
	int nmol = parameter.nmol;
	for (i = 0; i < nsfg; i++){
		symf_multi_vec_sa(molec, parameter, parameter.symftypeLst[i], i);
		for(int mi = 0; mi < nmol;mi++)
			atoms[molec[mi].center].sfg[i] = molec[mi].symf;
	}
	//calculate NN output and derivatives
	get_NN(atoms,parameter.center,natoms,NNvar);
}

//calculate symf and/or NN. Have option from INPUT file. stand alone only. Potentially faster than get_vec_symmetry_NN_vec_sa.  
void DAFED::get_vec_symmetry_NN_vec_sa_short_test(CAtom *atoms, Cmolpoint *molec,int &natoms, CParameter &parameter, NNvariables &NNvar, ofstream& symf_file)
{
	symf_multi_vec_sa_short_test(atoms,molec,parameter);
	//calculate NN output and derivatives
	
	if (parameter.onlysymcal)
		get_vec_symmetry_output(atoms, molec,natoms, parameter, symf_file);
	else
	{
		if(parameter.sym_new_flag)get_NN_sym_new(atoms,parameter.center,natoms,NNvar);
		else get_NN(atoms,parameter.center,natoms,NNvar);
	}
}

//printing out symf result to symf_result.txt
void DAFED::get_vec_symmetry_output(CAtom *atoms, Cmolpoint *molec,int &natoms, CParameter &parameter, ofstream& symf_file)
{
	unsigned short int mi, i;
	unsigned short int nsfg=parameter.nsfg;	
	unsigned short int nmol = parameter.nmol;
	for(mi=0;mi<nmol;++mi)       //loop over all molecules
	{
		for (i = 0; i < nsfg; ++i)
			symf_file << atoms[molec[mi].center].sfg[i] << " ";
		symf_file << endl;
	}

}



void DAFED::get_NNrdf(CAtom *atoms, Cmolpoint *molec,int &natoms, CParameter &parameter, double *box)
{
	unsigned short int mi;
	unsigned short int nmol = parameter.nmol;
	double tmpd, d;
	double diffx,diffy,diffz;
	for (mi = 0; mi < nmol; mi++)
	{
		d=20000000;
		if (molec[mi].n_neighbors < parameter.surface)		//if it is surface, 0
			d = 0; 
		else
		{
			for (int mj = 0; mj< nmol;mj++)
			{
				if (mi==mj)continue;	//skip if it is the same molecule
				if(!(molec[mj].n_neighbors < parameter.surface))continue;	//skip if they are not in surface
				tmpd = get_absDistance(molec[mi].center, molec[mj].center,diffx,diffy,diffz,atoms,box,parameter.triclinic);
				if(tmpd< d)
					d = tmpd;
			}
		}
		
		//This will check if I have molecule that has different number of atoms than natm. Also store the distance from surface to atom class.
		bool second_mol = false;
		if (molec[mi].atom_ID[parameter.natm-1]==-1)second_mol = true;
		for (int ti = 0; ti < parameter.natm;ti++)
		{
			if (molec[mi].atom_ID[ti] ==-1)break;
			atoms[molec[mi].atom_ID[ti]].d_surface = d;
			atoms[molec[mi].atom_ID[ti]].second_mol = second_mol;
		}
	}
}


//In this function, we calculate if the atom was classified as undefined or not wih NN output result. 
void DAFED::get_print_NNrdf(CAtom *atoms, Cmolpoint *molec,int &natoms, CParameter &parameter, NNvariables &NNvar, ofstream& file)
{
	unsigned short int mi;
        unsigned short int nmol = parameter.nmol;
	int flag;	//DK:chekcing if it is undefined or not 0: none of NNout are biggger than 0.5. 1: one of them are bigger than 0.5, 2: more than 2 are bigger than 0.5
	bool udf;	//This will be 0 or 1
	for (mi = 0; mi < nmol; mi++)
        {
		flag = 0;
		udf = 1;
		file << mi << " " << molec[mi].n_neighbors << " " << atoms[molec[mi].center].d_surface<< " ";
		for (int j =0;j < NNvar.output;j++)
		{
			file << atoms[molec[mi].center].NNout[j] << " ";
			if (flag == 0 && atoms[molec[mi].center].NNout[j] > 0.5)
			{
				flag = 1;
				udf = 0;
			}
			else if (flag == 1 && atoms[molec[mi].center].NNout[j] > 0.5)
			{
				flag =2;
				udf = 1;
			}
		}
		file << udf << " " << endl;
	}
}











//stand alone only. proto type of calculating neighbors. Still works.  
void DAFED::get_allNeighbourDistances_sa(CAtom *atoms, Cmolpoint *molec, int &natoms, double *box, CParameter &parameter)
{
	double d;		//absolute distance
	double diffx,diffy,diffz;	//xyz distance from target to neighbour 
	double fc, dfc;			//this is for cutoff and derivative of cutoff
	int nop = natoms;		//number of real atom in the proc
	int COvectype[2];			//data type of C and O
	int NNvectype[2];			//data type of N1 and N20	
	double rmin = parameter.rmin0;		//cutoff min
	double rmax = parameter.rmax0;		//cutoff max
	COvectype[0] = parameter.COvectype[0];	//assigning datatype 
	COvectype[1] = parameter.COvectype[1];	//assigning datatype
	NNvectype[0] = parameter.NNvectype[0];  //assigning datatype 
	NNvectype[1] = parameter.NNvectype[1];	//assigning datatype
	int nsfg=parameter.nsfg;		//assigning number of symmetry function
	int atm1_1 = -1, atm1_2=-1, atm2_1=-1, atm2_2=-1;	//initializing the atom id
	int molcount = 0;			//count number of mol in each
	int natm  = parameter.natm;
	
	////validate the size of the box
	if (rmax+2 > 0.5*box[0] || rmax+2 > 0.5*box[1] || rmax+2 > 0.5*box[2])			
	{
		cout << "cutoff(cutoff+skin(2)) is bigger than 1/2 of the box size. " << endl;
		cout << "box(x) = " << box[0] << endl;
		cout << "box(y) = " << box[1] << endl;
		cout << "box(z) = " << box[2] << endl;
		exit(1);
	}	
	//initializing atom class parameters. 
	for (int ti = 0;ti<nop;ti++)	//loop over all atoms
	{
		//DK: This is the part store all the symmetry function result that will be used for NN
		for (int k = 0;k<3;k++){	//loop over xyz
			atoms[ti].COvec[k] = 0.0;
			atoms[ti].NNvec[k] = 0.0;

		}			//loop over xyz
		for (int tsym = 0; tsym < nsfg; tsym++)	//loop over number of symmetry function
			atoms[ti].sfg[tsym] = 0;
		//assigning the center of the molecule and counting the number of mole. 
		if (atoms[ti].eletype == parameter.center)	
		{
			molec[molcount].center = ti;
			molcount++;
		}
	}				//loop over all atoms
	int atmcount[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	int non_imp_count[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	parameter.nmol = molcount;
	int dcount[molcount];	//using this for molecule count
	//initializing molecule class parameters
	for (int ii = 0;ii<molcount;ii++)	//loop over total mol in the proc
	{
		molec[ii].G2vecCO = 0;
		molec[ii].G2vecNN = 0;
		molec[ii].G3vecCO = 0;
		molec[ii].G3vecNN = 0;
		molec[ii].point2 = 0;
		molec[ii].point3 = 0;
		molec[ii].n_neighbors = 0;
		for (int k = 0;k<3;k++){	//loop over xyz for vec
			molec[ii].COvec[k] = 0;
			molec[ii].NNvec[k] = 0;
		}				//end of the loop over xyz for vec
		for (int tatm = 0;tatm<parameter.natm;tatm++)
			molec[ii].atom_ID[tatm] = -1;		//ID 0 --> C, 1--> O, 2--> N1, 3--> N2 4-7 -->H
		for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)	//loop for other parameters that need at least same number as 
		{
			molec[ii].fcv[tn] = 0.0;
			molec[ii].neighbors[tn] = -1;
			molec[ii].neighdist[tn] = 0;
			for (int k = 0;k<3;k++){	//loop over xyz for other parameters
				molec[ii].dfcvdx[tn][k] = 0.0;
				molec[ii].diff[tn][k] = 0;
			}//end loop for xyz
		}//end loop for other parameters
	}//end of total mol loop
	parameter.npairs = 0;
	
	for (int moli = 0; moli < molcount; moli++)
	{
		atmcount[moli] = 1;				//atmcount is 1D array and it should have 4 at the end. 
		non_imp_count[moli] = 4;
		dcount[moli]=0;
		int tii = molec[moli].center;			//get the center info, in the other words, this is getting the local id of carbon



		//I am checking if the center is one of the element for the vector or not. If not, center atom will be added to the end of the atom_ID list. 
		if(atoms[tii].eletype == COvectype[0])
			molec[moli].atom_ID[0] = tii;
		else if(atoms[tii].eletype == COvectype[1])
			molec[moli].atom_ID[1] = tii;
		else if(atoms[tii].eletype == NNvectype[0])
			molec[moli].atom_ID[2] = tii;
		else if(atoms[tii].eletype == NNvectype[1])
			molec[moli].atom_ID[3] = tii;
		else
		{
			molec[moli].atom_ID[natm-1] = tii; //if center atom isnot used for the vector calculation, ID will be added at the end of the atom_ID array
			non_imp_count[moli]++;
		}




		for (int tjj=0; tjj<nop; tjj++)		//loop over all neighbors
		{	
			if(tii==tjj)continue;			//if it is a same atom, skip it. 
			d = get_absDistance(tii,tjj,diffx,diffy,diffz,atoms,box,parameter.triclinic);
			if(d > rmax)continue;
			

			if (atoms[tjj].eletype == parameter.center) 		//pass only carbon/center of the mass here
			{
				if (d <= rmax){
					molec[moli].neighbors[dcount[moli]] = atoms[tjj].moltype;   //neighbor mol num
					molec[moli].neighdist[dcount[moli]] = d;
				
					molec[moli].diff[dcount[moli]][0] = diffx;
					molec[moli].diff[dcount[moli]][1] = diffy;
					molec[moli].diff[dcount[moli]][2] = diffz;
					if (abs(d) <= rmin){
						molec[moli].fcv[dcount[moli]] = 1.0;
					}
					else if (rmin < abs(d) && abs(d) <= rmax){					                                        
						fc = 0.5*(cos((d-rmin)*M_PI/(rmax-rmin))+1.0);
						molec[moli].fcv[dcount[moli]] = fc;                 //correct
						dfc = -0.5*M_PI/(rmax-rmin)*sin((d-rmin)*M_PI/(rmax-rmin));
						for (int k = 0;k<3;k++)
						{
							molec[moli].dfcvdx[dcount[moli]][k] = -dfc*molec[moli].diff[dcount[moli]][k]/d;	
						}
					}
					dcount[moli] +=1;
					
				}
			}
			
			if (atmcount[moli]==natm)continue;		//If you have all the atom info for this molecule, no need to search others
			if (atoms[tii].moltype != atoms[tjj].moltype)continue;          //If it is not in the same mol or C, I will skip
			if (atoms[tjj].eletype ==COvectype[0] )				//If this is C
			{
				molec[moli].atom_ID[0] = tjj;
				atmcount[moli]++;
			}
			if (atoms[tjj].eletype ==COvectype[1] )				//If this is O
			{
				molec[moli].atom_ID[1] = tjj;
				atmcount[moli]++;
			}
			else if (atoms[tjj].eletype ==NNvectype[0])			//If this is N1
			{
				molec[moli].atom_ID[2] = tjj;	
				atmcount[moli]++;
			}
			else if (atoms[tjj].eletype ==NNvectype[1])			//If this is N2
			{
				molec[moli].atom_ID[3] = tjj;	
				atmcount[moli]++;
			}
			else								//Technically, you dont need this info, but I keep this here for the future. 
			{
				molec[moli].atom_ID[non_imp_count[moli]++] = tjj;
				atmcount[moli]++;
			}
		}
	}
	//Validation and assiging CO and NN vector for the own molecule.  
	for (int mi = 0;mi<molcount;mi++){				//Validation for molec and getting atoms' molvec info.
		if (atmcount[mi] != natm)	
		{
			cout << "Warning: some of the atoms in the molecule are not detected!" << endl;	
			cout << "molec[mi].center: " << molec[mi].center << endl;
			cout << "atmcount[mi] : " << atmcount[mi] << endl;
			exit(1);
		}
		if (non_imp_count[mi] != natm)
		{
			cout << "Warning: some of the H/non_important atom in the molecule are not detected!" << endl;
			cout << "mi " << mi << endl;
			cout << "molec[mi].center: " << molec[mi].center << endl;
			cout << "non_imp_count[mi] : " << non_imp_count[mi] << endl;
			exit(1);
		}
		molec[mi].n_neighbors = dcount[mi];             //DK: this is for assignening number of neighbors
		if (dcount[mi]>MAXNUMBEROFNEIGHBORS)
		{
			cout << "dcount[mi] = "  << dcount[mi] << endl;
			cout << "You will have arrocation error. MAX you can hold is " << MAXNUMBEROFNEIGHBORS << " and now you have dcount[mi] " << dcount[mi] << endl;
			cout<<"Exiting programme..."<<endl;
			exit(1);
		}
		atm1_1 = molec[mi].atom_ID[0];
		atm1_2 = molec[mi].atom_ID[1];
		atm2_1 = molec[mi].atom_ID[2];
		atm2_2 = molec[mi].atom_ID[3];
		get_molvec_pbc(atm1_1,atm1_2,atoms[molec[mi].center].COvec, atoms,box,parameter.triclinic);
 		get_molvec_pbc(atm2_1,atm2_2,atoms[molec[mi].center].NNvec, atoms,box,parameter.triclinic);
		get_molvec_pbc(atm1_1,atm1_2,molec[mi].COvec, atoms,box,parameter.triclinic);
		get_molvec_pbc(atm2_1,atm2_2,molec[mi].NNvec, atoms,box,parameter.triclinic);
	}
}

//stand alone only. calculates the neighbor information here
void DAFED::get_allNeighbourDistances_gen_vec_sa(CAtom *atoms, Cmolpoint *molec, int &natoms, double *box, CParameter &parameter)
{
	double d;		//absolute distance
	double diffx,diffy,diffz;	//xyz distance from target to neighbour 
	double fc, dfc;			//this is for cutoff and derivative of cutoff
	int nop = natoms;		//number of real atom in the proc
	double rmin = parameter.rmin0;		//cutoff min
	double rmax = parameter.rmax0;		//cutoff max
	int nsfg=parameter.nsfg;		//assigning number of symmetry function
	int atm1 = -1, atm2=-1;	//initializing the atom id
	int molcount = 0;			//count number of mol in each
	int natm  = parameter.natm;
	
	////validate the size of the box
	if (rmax+2 > 0.5*box[0] || rmax+2 > 0.5*box[1] || rmax+2 > 0.5*box[2])			
	{
		cout << "cutoff(cutoff+skin(2)) is bigger than 1/2 of the box size. " << endl;
		cout << "box(x) = " << box[0] << endl;
		cout << "box(y) = " << box[1] << endl;
		cout << "box(z) = " << box[2] << endl;
		exit(1);
	}	
	//initializing atom class parameters. 
	for (int ti = 0;ti<nop;ti++)	//loop over all atoms
	{
		//DK: This is the part store all the symmetry function result that will be used for NN
		for (int k = 0;k<3;k++){	//loop over xyz
			atoms[ti].COvec[k] = 0.0;
			atoms[ti].NNvec[k] = 0.0;

		}			//loop over xyz
		for (int tsym = 0; tsym < nsfg; tsym++)	//loop over number of symmetry function
			atoms[ti].sfg[tsym] = 0;
		//assigning the center of the molecule and counting the number of mole. 
		if (atoms[ti].eletype == parameter.center)	
		{
			molec[molcount].center = ti;
			atoms[ti].sa_mol_id = molcount;
			molcount++;
		}
	}	
	//loop over all atoms
	int atmcount[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	parameter.nmol = molcount;
	int dcount[molcount];	//using this for molecule count
	//initializing molecule class parameters
	for (int ii = 0;ii<molcount;ii++)	//loop over total mol in the proc
	{
		molec[ii].G2vecCO = 0;
		molec[ii].G2vecNN = 0;
		molec[ii].G3vecCO = 0;
		molec[ii].G3vecNN = 0;
		molec[ii].point2 = 0;
		molec[ii].point3 = 0;
		molec[ii].n_neighbors = 0;
		for (int k = 0;k<3;k++){	//loop over xyz for vec
			molec[ii].COvec[k] = 0;
			molec[ii].NNvec[k] = 0;
		}				//end of the loop over xyz for vec
		for (int tatm = 0;tatm<parameter.natm;tatm++)
			molec[ii].atom_ID[tatm] = -1;		//ID 0 --> C, 1--> O, 2--> N1, 3--> N2 4-7 -->H
		for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)	//loop for other parameters that need at least same number as 
		{
			molec[ii].fcv[tn] = 0.0;
			molec[ii].neighbors[tn] = -1;
			molec[ii].neighdist[tn] = 0;
			for (int k = 0;k<3;k++){	//loop over xyz for other parameters
				molec[ii].dfcvdx[tn][k] = 0.0;
				molec[ii].diff[tn][k] = 0;
			}//end loop for xyz
		}//end loop for other parameters
	}//end of total mol loop
	parameter.npairs = 0;
	
	for (int moli = 0; moli < molcount; moli++)
	{
		atmcount[moli] = 1;				//atmcount is 1D array and it should have 4 at the end. 
		dcount[moli]=0;
		int tii = molec[moli].center;			//get the center info, in the other words, this is getting the local id of carbon

		//I am checking if the center is one of the element for the vector or not. If not, center atom will be added to the end of the atom_ID list. 
		molec[moli].atom_ID[0] = tii;
		molec[moli].vector_ID[0] = tii;			//point representation 
		molec[moli].vector_ID[1] = tii;			//this is dammy, but for the point representation.




		for (int tjj=0; tjj<nop; tjj++)		//loop over all neighbors
		{	
			if(tii==tjj)continue;			//if it is a same atom, skip it. 
			

			if (atoms[tjj].eletype == parameter.center) 		//pass only carbon/center of the mass here
			{
				d = get_absDistance(tii,tjj,diffx,diffy,diffz,atoms,box,parameter.triclinic);
				if(d > rmax)continue;
				if (d <= rmax){
					molec[moli].neighbors[dcount[moli]] = atoms[tjj].sa_mol_id;   //neighbor mol num
					molec[moli].neighdist[dcount[moli]] = d;
				
					molec[moli].diff[dcount[moli]][0] = diffx;
					molec[moli].diff[dcount[moli]][1] = diffy;
					molec[moli].diff[dcount[moli]][2] = diffz;
					if (abs(d) <= rmin){
						molec[moli].fcv[dcount[moli]] = 1.0;
					}
					else if (rmin < abs(d) && abs(d) <= rmax){					                                        
						fc = 0.5*(cos((d-rmin)*M_PI/(rmax-rmin))+1.0);
						molec[moli].fcv[dcount[moli]] = fc;                 //correct
						dfc = -0.5*M_PI/(rmax-rmin)*sin((d-rmin)*M_PI/(rmax-rmin));
						for (int k = 0;k<3;k++)
						{
							molec[moli].dfcvdx[dcount[moli]][k] = -dfc*molec[moli].diff[dcount[moli]][k]/d;	
						}
					}
					dcount[moli] +=1;
					
				}
			}
			
			if (atmcount[moli]==natm)continue;		//If you have all the atom info for this molecule, no need to search others
			if (atoms[tii].moltype != atoms[tjj].moltype)continue;          //If it is not in the same mol or C, I will skip
			
			molec[moli].atom_ID[atmcount[moli]] = tjj;
			atmcount[moli]++;
		}
		for (int sfgt = 2; sfgt < 2*parameter.symftype; sfgt++)
		{
			for (int ti = 0; ti < natm; ti++)
			{
				if(molec[moli].atom_ID[ti]==-1)break;                   //If it does not exist, we skip it. 
				if (parameter.vectype[sfgt] == atoms[molec[moli].atom_ID[ti]].eletype)
				{
					molec[moli].vector_ID[sfgt] = molec[moli].atom_ID[ti];
					break;
				}
			}
		}
	}
	//Validation and assiging CO and NN vector for the own molecule.  
	for (int mi = 0;mi<molcount;mi++){				//Validation for molec and getting atoms' molvec info.
		//this used to be a fail safe program for bulk system. If you have defect, this part need to be commented out, otherwise, put thes codes back back  
		//if (atmcount[mi] != natm)	
		//{
		//	cout << "Warning: some of the atoms in the molecule are not detected!" << endl;	
		//	cout << "molec[mi].center: " << molec[mi].center << endl;
		//	cout << "atmcount[mi] : " << atmcount[mi] << endl;
		//	exit(1);
		//}
		molec[mi].n_neighbors = dcount[mi];             //DK: this is for assignening number of neighbors
		for (int ti = 0; ti < natm;ti++)
			if(molec[mi].atom_ID[ti] !=-1)		//DK: due to defects, some of atoms will be missing. I need to take care of that.
				atoms[molec[mi].atom_ID[ti]].n_neighbors = dcount[mi];
		if (dcount[mi]>MAXNUMBEROFNEIGHBORS)
		{
			cout << "dcount[mi] = "  << dcount[mi] << endl;
			cout << "You will have arrocation error. MAX you can hold is " << MAXNUMBEROFNEIGHBORS << " and now you have dcount[mi] " << dcount[mi] << endl;
			cout<<"Exiting programme..."<<endl;
			exit(1);
		}

		for (int sfg = 1; sfg < parameter.symftype; sfg++)
		{
			atm1 = molec[mi].vector_ID[sfg*2];
			atm2 = molec[mi].vector_ID[sfg*2+1];
			get_molvec_pbc(atm1,atm2,molec[mi].vecpoint[sfg], atoms,box,parameter.triclinic);
		}
		
	}
}


//stand alone only. calculates the neighbor information here
void DAFED::get_allNeighbourDistances_sa_sym_new(CAtom *atoms, Cmolpoint *molec, int &natoms, double *box, CParameter &parameter)
{
	double d;		//absolute distance
	double diffx,diffy,diffz;	//xyz distance from target to neighbour 
	double fc, dfc;			//this is for cutoff and derivative of cutoff
	int nop = natoms;		//number of real atom in the proc
	double rmin = parameter.rmin0;		//cutoff min
	double rmax = parameter.rmax0;		//cutoff max
	int nsfg=parameter.nsfg;		//assigning number of symmetry function
	int atm1 = -1, atm2=-1;	//initializing the atom id
	int molcount = 0;			//count number of mol in each
	int natm  = parameter.natm;
	
	////validate the size of the box
	if (rmax+2 > 0.5*box[0] || rmax+2 > 0.5*box[1] || rmax+2 > 0.5*box[2])			
	{
		cout << "cutoff(cutoff+skin(2)) is bigger than 1/2 of the box size. " << endl;
		cout << "box(x) = " << box[0] << endl;
		cout << "box(y) = " << box[1] << endl;
		cout << "box(z) = " << box[2] << endl;
		exit(1);
	}	
	//initializing atom class parameters. 
	for (int ti = 0;ti<nop;ti++)	//loop over all atoms
	{
		//DK: This is the part store all the symmetry function result that will be used for NN
		for (int k = 0;k<3;k++){	//loop over xyz
			atoms[ti].COvec[k] = 0.0;
			atoms[ti].NNvec[k] = 0.0;

		}			//loop over xyz
		for (int tsym = 0; tsym < nsfg; tsym++)	//loop over number of symmetry function
			atoms[ti].sfg[tsym] = 0;
		//assigning the center of the molecule and counting the number of mole. 
		if (atoms[ti].symtype == parameter.center)	
		{
			molec[molcount].center = ti;
			atoms[ti].sa_mol_id = molcount;
			molcount++;
		}
	}	
	//loop over all atoms
	int atmcount[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	parameter.nmol = molcount;
	int dcount[molcount];	//using this for molecule count
	//initializing molecule class parameters
	for (int ii = 0;ii<molcount;ii++)	//loop over total mol in the proc
	{
		molec[ii].G2vecCO = 0;
		molec[ii].G2vecNN = 0;
		molec[ii].G3vecCO = 0;
		molec[ii].G3vecNN = 0;
		molec[ii].point2 = 0;
		molec[ii].point3 = 0;
		molec[ii].n_neighbors = 0;
		for (int k = 0;k<3;k++){	//loop over xyz for vec
			molec[ii].COvec[k] = 0;
			molec[ii].NNvec[k] = 0;
		}				//end of the loop over xyz for vec
		for (int tatm = 0;tatm<parameter.natm;tatm++)
			molec[ii].atom_ID[tatm] = -1;		//ID 0 --> C, 1--> O, 2--> N1, 3--> N2 4-7 -->H
		for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)	//loop for other parameters that need at least same number as 
		{
			molec[ii].fcv[tn] = 0.0;
			molec[ii].neighbors[tn] = -1;
			molec[ii].neighdist[tn] = 0;
			for (int k = 0;k<3;k++){	//loop over xyz for other parameters
				molec[ii].dfcvdx[tn][k] = 0.0;
				molec[ii].diff[tn][k] = 0;
			}//end loop for xyz
		}//end loop for other parameters
	}//end of total mol loop
	parameter.npairs = 0;
	
	for (int moli = 0; moli < molcount; moli++)
	{
		atmcount[moli] = 1;				//atmcount is 1D array and it should have 4 at the end. 
		dcount[moli]=0;
		int tii = molec[moli].center;			//get the center info, in the other words, this is getting the local id of carbon

		//I am checking if the center is one of the element for the vector or not. If not, center atom will be added to the end of the atom_ID list. 
		molec[moli].atom_ID[0] = tii;
		molec[moli].vector_ID[0] = tii;			//point representation 
		molec[moli].vector_ID[1] = tii;			//this is dammy, but for the point representation.




		for (int tjj=0; tjj<nop; tjj++)		//loop over all neighbors
		{	
			if(tii==tjj)continue;			//if it is a same atom, skip it. 
			

			if (atoms[tjj].symtype == parameter.center) 		//pass only center atom here
			{
				d = get_absDistance(tii,tjj,diffx,diffy,diffz,atoms,box,parameter.triclinic);
				if(d > rmax)continue;
				if (d <= rmax){
					molec[moli].neighbors[dcount[moli]] = atoms[tjj].sa_mol_id;   //neighbor mol num
					molec[moli].neighdist[dcount[moli]] = d;
				
					molec[moli].diff[dcount[moli]][0] = diffx;
					molec[moli].diff[dcount[moli]][1] = diffy;
					molec[moli].diff[dcount[moli]][2] = diffz;
					if (abs(d) <= rmin){
						molec[moli].fcv[dcount[moli]] = 1.0;
					}
					else if (rmin < abs(d) && abs(d) <= rmax){					                                        
						fc = 0.5*(cos((d-rmin)*M_PI/(rmax-rmin))+1.0);
						molec[moli].fcv[dcount[moli]] = fc;                 //correct
						dfc = -0.5*M_PI/(rmax-rmin)*sin((d-rmin)*M_PI/(rmax-rmin));
						for (int k = 0;k<3;k++)
						{
							molec[moli].dfcvdx[dcount[moli]][k] = -dfc*molec[moli].diff[dcount[moli]][k]/d;	
						}
					}
					dcount[moli] +=1;
					
				}
			}
			
			if (atmcount[moli]==natm)continue;		//If you have all the atom info for this molecule, no need to search others
			if (atoms[tii].moltype != atoms[tjj].moltype)continue;          //If it is not in the same mol or C, I will skip
			
			molec[moli].atom_ID[atmcount[moli]] = tjj;
			atmcount[moli]++;
		}
		for (int sfgt = 2; sfgt < 2*parameter.symftype; sfgt++)
		{
			for (int ti = 0; ti < natm; ti++)
			{
				if(molec[moli].atom_ID[ti]==-1)break;                   //If it does not exist, we skip it. 
				if (parameter.vectype[sfgt] == atoms[molec[moli].atom_ID[ti]].symtype)	//we assign the vector here
				{
					molec[moli].vector_ID[sfgt] = molec[moli].atom_ID[ti];
					break;
				}
			}
		}
	}
	//Validation and assiging CO and NN vector for the own molecule.  
	for (int mi = 0;mi<molcount;mi++){				//Validation for molec and getting atoms' molvec info.
		//this used to be a fail safe program for bulk system. If you have defect, this part need to be commented out, otherwise, put thes codes back back  
		//if (atmcount[mi] != natm)	
		//{
		//	cout << "Warning: some of the atoms in the molecule are not detected!" << endl;	
		//	cout << "molec[mi].center: " << molec[mi].center << endl;
		//	cout << "atmcount[mi] : " << atmcount[mi] << endl;
		//	exit(1);
		//}
		molec[mi].n_neighbors = dcount[mi];             //DK: this is for assignening number of neighbors
		for (int ti = 0; ti < natm;ti++)
			if(molec[mi].atom_ID[ti] !=-1)		//DK: due to defects, some of atoms will be missing. I need to take care of that.
				atoms[molec[mi].atom_ID[ti]].n_neighbors = dcount[mi];
		if (dcount[mi]>MAXNUMBEROFNEIGHBORS)
		{
			cout << "dcount[mi] = "  << dcount[mi] << endl;
			cout << "You will have arrocation error. MAX you can hold is " << MAXNUMBEROFNEIGHBORS << " and now you have dcount[mi] " << dcount[mi] << endl;
			cout<<"Exiting programme..."<<endl;
			exit(1);
		}

		for (int sfg = 1; sfg < parameter.symftype; sfg++)
		{
			atm1 = molec[mi].vector_ID[sfg*2];
			atm2 = molec[mi].vector_ID[sfg*2+1];
			get_molvec_pbc(atm1,atm2,molec[mi].vecpoint[sfg], atoms,box,parameter.triclinic);
		}
		
	}
}




//lammps only. proto type. Can hold only two vectors
//This is making neighbours information. 
//We use the lammps neighbor information. 
void DAFED::get_allNeighbourDistances_lmp_new(CAtom *atoms, Cmolpoint *molec, CNeighvariables &lmpneigh, int &natoms, double *box, CParameter &parameter)
{
	double d;		//absolute distance
	double diffx,diffy,diffz;	//xyz distance from target to neighbour 
	double fc, dfc;			//this is for cutoff and derivative of cutoff
	int jnum;			//number of neighbours for the specific target. 
	int *jlist;			//neighbour list
	int nop = lmpneigh.inum;		//number of real atom in the proc
	int COvectype[2];			//data type of C and O
	int NNvectype[2];			//data type of N1 and N20	
	double rmin = parameter.rmin0;		//cutoff min
	double rmax = parameter.rmax0;		//cutoff max
	COvectype[0] = parameter.COvectype[0];	//assigning datatype 
	COvectype[1] = parameter.COvectype[1];	//assigning datatype
	NNvectype[0] = parameter.NNvectype[0];  //assigning datatype 
	NNvectype[1] = parameter.NNvectype[1];	//assigning datatype
	int nsfg=parameter.nsfg;		//assigning number of symmetry function
	int atmC = -1, atmO=-1, atmN1=-1, atmN2=-1;	//initializing the atom id
	int molcount = 0;			//count number of mol in each

	
	////validate the size of the box
	if (rmax+2 > 0.5*box[0] || rmax+2 > 0.5*box[1] || rmax+2 > 0.5*box[2])			
	{
		cout << "cutoff(cutoff+skin(2)) is bigger than 1/2 of the box size. " << endl;
		cout << "box(x) = " << box[0] << endl;
		cout << "box(y) = " << box[1] << endl;
		cout << "box(z) = " << box[2] << endl;
		exit(1);
	}	
	
	//initializing atom class parameters. 
	for (int ti = 0;ti<nop;ti++)	//loop over all atoms
	{
		//DK: This is the part store all the symmetry function result that will be used for NN
		for (int k = 0;k<3;k++){	//loop over xyz
			atoms[ti].COvec[k] = 0.0;
			atoms[ti].NNvec[k] = 0.0;

		}			//loop over xyz
		for (int tsym = 0; tsym < nsfg; tsym++)	//loop over number of symmetry function
			atoms[ti].sfg[tsym] = 0;
		//assigning the center of the molecule and counting the number of mole. 
		if (atoms[ti].eletype == parameter.center)	
		{	
			molec[molcount].center = ti;
			molcount++;
		}
	}				//loop over all atoms
	int atmcount[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	int Hcount[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	int lmpmol = molcount;	//number of molecule in this proc
	lmpneigh.lmpmol = lmpmol;	//number of mol in the proc need to be updated.

	int dcount[lmpmol];	//using this for molecule count
	//initializing molecule class parameters
	for (int ii = 0;ii<lmpmol;ii++)	//loop over total mol in the proc
	{
		molec[ii].G2vecCO = 0;
		molec[ii].G2vecNN = 0;
		molec[ii].G3vecCO = 0;
		molec[ii].G3vecNN = 0;
		molec[ii].point2 = 0;
		molec[ii].point3 = 0;
		molec[ii].n_neighbors = 0;
		for (int k = 0;k<3;k++){	//loop over xyz for vec
			molec[ii].COvec[k] = 0;
			molec[ii].NNvec[k] = 0;
		}				//end of the loop over xyz for vec
		for (int tatm = 0;tatm<parameter.natm;tatm++)
			molec[ii].atom_ID[tatm] = -1;		//ID 0 --> C, 1--> O, 2--> N1, 3--> N2 4-7 -->H
		for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)	//loop for other parameters that need at least same number as 
		{
			molec[ii].fcv[tn] = 0.0;
			molec[ii].neighbors[tn] = -1;
			molec[ii].neighdist[tn] = 0;
			for (int k = 0;k<3;k++){	//loop over xyz for other parameters
				molec[ii].dfcvdx[tn][k] = 0.0;
				molec[ii].diff[tn][k] = 0;
			}//end loop for xyz
		}//end loop for other parameters
	}//end of total mol loop
	parameter.npairs = 0;
	
	for (int moli = 0; moli < lmpmol; moli++)
	{
		atmcount[moli] = 1;				//atmcount is 1D array and it should have 4 at the end. 
		Hcount[moli] = 4;
		dcount[moli]=0;
		int tii = molec[moli].center;			//get the center info, in the other words, this is getting the local id of carbon
		int ti = lmpneigh.ilist[tii];			//get local id?
		jlist = lmpneigh.firstneigh[ti];		//list of ti neighbors 
		jnum = lmpneigh.numneigh[ti];			//number of ti neighbors for ti
		
		
		
		if(atoms[tii].eletype == COvectype[0])
			molec[moli].atom_ID[0] = tii;
		else if(atoms[tii].eletype == COvectype[1])
			molec[moli].atom_ID[1] = tii;
		else if(atoms[tii].eletype == NNvectype[0])
			molec[moli].atom_ID[2] = tii;
		else if(atoms[tii].eletype == NNvectype[1])
			molec[moli].atom_ID[3] = tii;
		else
		{
			molec[moli].atom_ID[parameter.natm-1] = tii; //if center atom isnot used for the vector calculation, ID will be added at the end of the atom_ID array
			Hcount[moli]++;
		}
		
		
		for (int tjj=0; tjj<jnum; tjj++)		//loop over all neighbors
		{	
			int tj = jlist[tjj];	
			tj &= MY_NEIGHMASK;
			if (atoms[tj].eletype == parameter.center) 		//pass only carbon/center of the mass here
			{
				if (atoms[tj].moltype!=atoms[ti].moltype)
				{
					d = get_absDistance(ti,tj,diffx,diffy,diffz,atoms,box,parameter.triclinic);
					if (d <= rmax){
						molec[moli].neighbors[dcount[moli]] = tj;   //neighbor num
						molec[moli].neighdist[dcount[moli]] = d;
				
						molec[moli].diff[dcount[moli]][0] = diffx;
						molec[moli].diff[dcount[moli]][1] = diffy;
						molec[moli].diff[dcount[moli]][2] = diffz;
						if (abs(d) <= rmin){
							molec[moli].fcv[dcount[moli]] = 1.0;
						}
						else if (rmin < abs(d) && abs(d) <= rmax){					                                        
							fc = 0.5*(cos((d-rmin)*M_PI/(rmax-rmin))+1.0);
							molec[moli].fcv[dcount[moli]] = fc;                 //correct
							dfc = -0.5*M_PI/(rmax-rmin)*sin((d-rmin)*M_PI/(rmax-rmin));
							for (int k = 0;k<3;k++)
							{
								molec[moli].dfcvdx[dcount[moli]][k] = -dfc*molec[moli].diff[dcount[moli]][k]/d;	
							}
						}
						dcount[moli] +=1;
					}
				}
			}
			
			if (atmcount[moli]==parameter.natm)continue;		//If you have all the atom info for this molecule, no need to search others
			if (atoms[ti].moltype != atoms[tj].moltype || atoms[tj].eletype == parameter.center)continue;          //If it is not in the same mol or C, I will skip
			if(atoms[ti].moltype == 0){
			}
			if (atoms[tj].eletype ==COvectype[0] )				//If this is C
			{
				molec[moli].atom_ID[0] = tj;
				atmcount[moli]++;
			}
			else if (atoms[tj].eletype ==COvectype[1] )				//If this is O
			{
				molec[moli].atom_ID[1] = tj;
				atmcount[moli]++;
			}
			else if (atoms[tj].eletype ==NNvectype[0])			//If this is N1
			{
				molec[moli].atom_ID[2] = tj;	
				atmcount[moli]++;
			}
			else if (atoms[tj].eletype ==NNvectype[1])			//If this is N2
			{
				molec[moli].atom_ID[3] = tj;	
				atmcount[moli]++;
			}
			else								//Technically, you dont need this info, but I keep this here for the future. 
			{
				molec[moli].atom_ID[Hcount[moli]++] = tj;
				atmcount[moli]++;
			}
		}
	}
	//Validation and assiging CO and NN vector for the own molecule.  
	for (int mi = 0;mi<lmpmol;mi++){				//Validation for molec and getting atoms' molvec info.
		
		if (atmcount[mi] != parameter.natm)	
		{
			cout << "Warning: some of the atoms in the molecule are not detected!" << endl;	
			cout << "molec[mili].center: " << molec[mi].center << endl;
			cout << "atmcount[mi] : " << atmcount[mi] << endl;
			for (int at = 0; at<parameter.natm; at++)
				cout << molec[mi].atom_ID[at] << endl;
			exit(1);
		}
		if (Hcount[mi] != parameter.natm)
		{
			cout << "Warning: some of the H in the molecule are not detected!" << endl;
			cout << "molec_num: " << mi << endl;
			cout << "molec[mili].center: " << molec[mi].center << endl;
			cout << "atmcount[mi] : " << atmcount[mi] << endl;
			cout << "Hcount[mi] : " << Hcount[mi] << endl;
			for (int at = 0; at<parameter.natm; at++)
			{
				cout << "molec[mi].atom_ID[at]  = " << molec[mi].atom_ID[at] << endl;
				cout << "atoms[molec[mi].atom_ID[at]].eletype = " << atoms[molec[mi].atom_ID[at]].eletype << endl;
				cout << "atoms[molec[mi].atom_ID[at]].moltype = " << atoms[molec[mi].atom_ID[at]].moltype << endl;
				cout << "atoms[molec[mi].atom_ID[at]].pos[x] = " << atoms[molec[mi].atom_ID[at]].pos[0] << endl;
				cout << "atoms[molec[mi].atom_ID[at]].pos[y] = " << atoms[molec[mi].atom_ID[at]].pos[1] << endl;
				cout << "atoms[molec[mi].atom_ID[at]].pos[z] = " << atoms[molec[mi].atom_ID[at]].pos[2] << endl;
			}
			exit(1);
		}
		molec[mi].n_neighbors = dcount[mi];             //DK: this is for assignening number of neighbors
		if (dcount[mi]>MAXNUMBEROFNEIGHBORS)
		{
			cout << "dcount[mi] = "  << dcount[mi] << endl;
			cout << "You will have arrocation error. MAX you can hold is " << MAXNUMBEROFNEIGHBORS << " and now you have dcount[mi] " << dcount[mi] << endl;
			cout<<"Exiting programme..."<<endl;
			exit(1);
		}
		atmC = molec[mi].atom_ID[0];
		atmO = molec[mi].atom_ID[1];
		atmN1 = molec[mi].atom_ID[2];
		atmN2 = molec[mi].atom_ID[3];
		get_molvec_pbc(atmC,atmO,atoms[molec[mi].center].COvec, atoms,box,parameter.triclinic);
		get_molvec_pbc(atmN1,atmN2,atoms[molec[mi].center].NNvec, atoms,box,parameter.triclinic);
	}
}

//lammps only. calculate neighbor information. 
void DAFED::get_allNeighbourDistances_gen_vec_lmp(CAtom *atoms, Cmolpoint *molec, CNeighvariables &lmpneigh, int &natoms, double *box, CParameter &parameter)
{
	double d;		//absolute distance
	double diffx,diffy,diffz;	//xyz distance from target to neighbour 
	double fc, dfc;			//this is for cutoff and derivative of cutoff
        int jnum;                       //number of neighbours for the specific target. 
        int *jlist;                     //neighbour list
	int nop = lmpneigh.inum;		//number of real atom in the proc
	double rmin = parameter.rmin0;		//cutoff min
	double rmax = parameter.rmax0;		//cutoff max
	int nsfg=parameter.nsfg;		//assigning number of symmetry function
	int atm1 = -1, atm2=-1;	//initializing the atom id
	int molcount = 0;			//count number of mol in each
	int natm  = parameter.natm;
	
	////validate the size of the box
	if (rmax+2 > 0.5*box[0] || rmax+2 > 0.5*box[1] || rmax+2 > 0.5*box[2])			
	{
		cout << "cutoff(cutoff+skin(2)) is bigger than 1/2 of the box size. " << endl;
		cout << "box(x) = " << box[0] << endl;
		cout << "box(y) = " << box[1] << endl;
		cout << "box(z) = " << box[2] << endl;
		exit(1);
	}	
	//initializing atom class parameters. 
	for (int ti = 0;ti<nop;ti++)	//loop over all atoms
	{
		//DK: This is the part store all the symmetry function result that will be used for NN
		for (int k = 0;k<3;k++){	//loop over xyz
			atoms[ti].COvec[k] = 0.0;
			atoms[ti].NNvec[k] = 0.0;

		}			//loop over xyz
		for (int tsym = 0; tsym < nsfg; tsym++)	//loop over number of symmetry function
			atoms[ti].sfg[tsym] = 0;
		//assigning the center of the molecule and counting the number of mole. 
		if (atoms[ti].eletype == parameter.center)	
		{
			molec[molcount].center = ti;
			molcount++;
		}
	}				//loop over all atoms
	int atmcount[molcount];	//count number of important atoms in the molecule counted. For urea, it have to be 4 at the end
	lmpneigh.lmpmol = molcount;
	int dcount[molcount];	//using this for molecule count
	//initializing molecule class parameters
	for (int ii = 0;ii<molcount;ii++)	//loop over total mol in the proc
	{
		molec[ii].G2vecCO = 0;
		molec[ii].G2vecNN = 0;
		molec[ii].G3vecCO = 0;
		molec[ii].G3vecNN = 0;
		molec[ii].point2 = 0;
		molec[ii].point3 = 0;
		molec[ii].n_neighbors = 0;
		for (int tatm = 0;tatm<parameter.natm;tatm++)
			molec[ii].atom_ID[tatm] = -1;		//ID 0 --> C, 1--> O, 2--> N1, 3--> N2 4-7 -->H
		for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)	//loop for other parameters that need at least same number as 
		{
			molec[ii].fcv[tn] = 0.0;
			molec[ii].neighbors[tn] = -1;
			molec[ii].neighdist[tn] = 0;
			for (int k = 0;k<3;k++){	//loop over xyz for other parameters
				molec[ii].dfcvdx[tn][k] = 0.0;
				molec[ii].diff[tn][k] = 0;
			}//end loop for xyz
		}//end loop for other parameters
	}//end of total mol loop
	parameter.npairs = 0;
	
	for (int moli = 0; moli < molcount; moli++)
	{
		atmcount[moli] = 1;				//atmcount is 1D array and it should have 4 at the end. 
		dcount[moli]=0;
		int tii = molec[moli].center;			//get the center info, in the other words, this is getting the local id of carbon
		int ti = lmpneigh.ilist[tii];			//get local id? I think this is same as tii. Need to comfirm. 
		jlist = lmpneigh.firstneigh[ti];		//list of ti neighbors
		jnum = lmpneigh.numneigh[ti];                   //number of ti neighbors for ti


		//I am checking if the center is one of the element for the vector or not. If not, center atom will be added to the end of the atom_ID list. 
		molec[moli].atom_ID[0] = tii;
		molec[moli].vector_ID[0] = tii;			//point representation 
		molec[moli].vector_ID[1] = tii;			//this is dammy, but for the point representation.




		for (int tjj=0; tjj<jnum; tjj++)		//loop over all neighbors
		{	
			int tj = jlist[tjj];
                        tj &= MY_NEIGHMASK;
			if(tii==tj)continue;			//if it is a same atom, skip it. 
			

			if (atoms[tj].eletype == parameter.center) 		//pass only carbon/center of the mass here
			{
				d = get_absDistance(tii,tj,diffx,diffy,diffz,atoms,box,parameter.triclinic);
				if(d > rmax)continue;
				if (d <= rmax){
					molec[moli].neighbors[dcount[moli]] = tj;   //neighbour id
					molec[moli].neighdist[dcount[moli]] = d;
				
					molec[moli].diff[dcount[moli]][0] = diffx;
					molec[moli].diff[dcount[moli]][1] = diffy;
					molec[moli].diff[dcount[moli]][2] = diffz;
					if (abs(d) <= rmin){
						molec[moli].fcv[dcount[moli]] = 1.0;
					}
					else if (rmin < abs(d) && abs(d) <= rmax){					                                        
						fc = 0.5*(cos((d-rmin)*M_PI/(rmax-rmin))+1.0);
						molec[moli].fcv[dcount[moli]] = fc;                 //correct
						dfc = -0.5*M_PI/(rmax-rmin)*sin((d-rmin)*M_PI/(rmax-rmin));
						for (int k = 0;k<3;k++)
						{
							molec[moli].dfcvdx[dcount[moli]][k] = -dfc*molec[moli].diff[dcount[moli]][k]/d;	
						}
					}
					dcount[moli] +=1;
					
				}
			}
			
			if (atmcount[moli]==natm)continue;		//If you have all the atom info for this molecule, no need to search others
			if (atoms[tii].moltype != atoms[tj].moltype)continue;          //If it is not in the same mol or C, I will skip
			
		
			
			molec[moli].atom_ID[atmcount[moli]] = tj;
			atmcount[moli]++;
		}
		for (int sfgt = 2; sfgt < 2*parameter.symftype; sfgt++)
		{
			for (int ti = 0; ti < natm; ti++)
			{
				if(molec[moli].atom_ID[ti]==-1)break;			//If it does not exist, we skip it. 
				if (parameter.vectype[sfgt] == atoms[molec[moli].atom_ID[ti]].eletype)
				{
					molec[moli].vector_ID[sfgt] = molec[moli].atom_ID[ti];
					break;
				}
			}
		}
	}
	//Validation and assiging CO and NN vector for the own molecule.  
	for (int mi = 0;mi<molcount;mi++){				//Validation for molec and getting atoms' molvec info.
/*commented out because if you have defect that has different atomic number, it would not work. Otherwise, you can put these comments back. 
 		if (atmcount[mi] != natm)	
		{
			cout << "Warning: some of the atoms in the molecule are not detected!" << endl;	
			cout << "molec[mi].center: " << molec[mi].center << endl;
			cout << "atmcount[mi] : " << atmcount[mi] << endl;
			exit(1);
		}
*/		molec[mi].n_neighbors = dcount[mi];             //DK: this is for assignening number of neighbors
		if (dcount[mi]>MAXNUMBEROFNEIGHBORS)
		{
			cout << "dcount[mi] = "  << dcount[mi] << endl;
			cout << "You will have arrocation error. MAX you can hold is " << MAXNUMBEROFNEIGHBORS << " and now you have dcount[mi] " << dcount[mi] << endl;
			cout<<"Exiting programme..."<<endl;
			exit(1);
		}

		for (int sfg = 1; sfg < parameter.symftype; sfg++)
		{
			atm1 = molec[mi].vector_ID[sfg*2];
			atm2 = molec[mi].vector_ID[sfg*2+1];
			get_molvec_pbc(atm1,atm2,molec[mi].vecpoint[sfg], atoms,box,parameter.triclinic);
			get_molvec_pbc(atm1,atm2,atoms[molec[mi].center].vecpoint[sfg], atoms,box,parameter.triclinic);
		}
		
	}
}

//get distance between two atoms
double DAFED::get_absDistance(int ti ,int tj,double &diffx ,double &diffy,double &diffz, CAtom *molecules, double *box, bool triclinic)
{
	double abs;
	diffx = molecules[tj].pos[0] - molecules[ti].pos[0];
  	diffy = molecules[tj].pos[1] - molecules[ti].pos[1];
  	diffz = molecules[tj].pos[2] - molecules[ti].pos[2];
 	if (!(triclinic))
	{
  		//nearest image
  		if (diffx >  box[0]/2.0) {diffx = diffx - box[0];};
  		if (diffx < -box[0]/2.0) {diffx = diffx + box[0];};
  		if (diffy >  box[1]/2.0) {diffy = diffy - box[1];};
  		if (diffy < -box[1]/2.0) {diffy = diffy + box[1];};
  		if (diffz >  box[2]/2.0) {diffz = diffz - box[2];};
  		if (diffz < -box[2]/2.0) {diffz = diffz + box[2];};
	}

	else
	{
		//for z
		if (diffz > box[2]/2.0)
		{
			diffz -= box[2];	//diffz - zbox
			diffy -= box[3];	//diffy - yz
			diffx -= box[4];	//diffx - xz
		}
		if (diffz < -box[2]/2.0)
		{
			diffz += box[2];	//diffz - zbox
			diffy += box[3];	//diffy - yz
			diffx += box[4];	//diffx - xz
		}
		//for y
		if (diffy > box[1]/2.0)
		{
			diffy -= box[1];	//diffy - ybox
			diffx -= box[5];	//diffx - xy
		}
		if (diffy < -box[1]/2.0)
		{
			diffy += box[1];	//diffy - ybox
			diffx += box[5];	//diffx - xy
		}
		//for x
		if (diffx >  box[0]/2.0) {diffx = diffx - box[0];};
        	if (diffx < -box[0]/2.0) {diffx = diffx + box[0];};
	}
  	abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
  	return abs;
}

//moving ghost atoms 
void DAFED::get_molvec_pbc(int ti ,int tj,double *molvec, CAtom *atoms, double *box, bool triclinic)
{
	if (!(triclinic))
	{
		molvec[0] = atoms[tj].pos[0] - atoms[ti].pos[0];
		molvec[1] = atoms[tj].pos[1] - atoms[ti].pos[1];
		molvec[2] = atoms[tj].pos[2] - atoms[ti].pos[2];

  		if (molvec[0] >  box[0]/2.0) {molvec[0] = molvec[0] - box[0];};
  		if (molvec[0] < -box[0]/2.0) {molvec[0] = molvec[0] + box[0];};
  		if (molvec[1] >  box[1]/2.0) {molvec[1] = molvec[1] - box[1];};
  		if (molvec[1] < -box[1]/2.0) {molvec[1] = molvec[1] + box[1];};
  		if (molvec[2] >  box[2]/2.0) {molvec[2] = molvec[2] - box[2];};
  		if (molvec[2] < -box[2]/2.0) {molvec[2] = molvec[2] + box[2];};			
	}
	else
	{
		molvec[0] = atoms[tj].pos[0] - atoms[ti].pos[0];
		molvec[1] = atoms[tj].pos[1] - atoms[ti].pos[1];
		molvec[2] = atoms[tj].pos[2] - atoms[ti].pos[2];
		if (molvec[2] > box[2]/2.0)
		{
			//diffz = diffz - box[2] - box[3] - box[4]; 	//diffz - zbox - yz - xz	
			molvec[2] -= box[2];	//diffz - zbox
			molvec[1] -= box[3];	//diffy - yz
			molvec[0] -= box[4];	//diffx - xz
		}
		if (molvec[2] < -box[2]/2.0)
		{
			//diffz = diffz + box[2] + box[3] + box[4]; 	//diffz + zbox + yz + xz	
			molvec[2] += box[2];	//diffz - zbox
			molvec[1] += box[3];	//diffy - yz
			molvec[0] += box[4];	//diffx - xz
		}
		if (molvec[1] > box[1]/2.0)
		{
			//diffy = diffy - box[1] - box[5];	//diffy - ybox - xy	
			molvec[1] -= box[1];	//diffy - ybox
			molvec[0] -= box[5];	//diffx - xy
		}
		if (molvec[1] < -box[1]/2.0)
		{
			//diffy = diffy + box[1] + box[5]; 	//diffy + ybox + xy	
			molvec[1] += box[1];	//diffy - ybox
			molvec[0] += box[5];	//diffx - xy
		}
		if (molvec[0] >  box[0]/2.0) {molvec[0] = molvec[0] - box[0];};
        	if (molvec[0] < -box[0]/2.0) {molvec[0] = molvec[0] + box[0];};	
	}
}

//------------------------------Class CMolecules---------------------------------
//Constructor for atoms arrays
DAFED::CAtom::CAtom()
{
}


//Destructor
DAFED::CAtom::~CAtom()
{
}

DAFED::Cmolpoint::Cmolpoint()
{
}


DAFED::Cmolpoint::~Cmolpoint()
{
}

//----------------------------Class CParameter------------------------------------

//Constructor
DAFED::CParameter::CParameter()
{
	//this->neighbourdistance = -1.0;
	this->rmin0 = -1.0;
	this->rmax0 = -1.0;
	this->rmin1 = -1.0;
	this->rmax1 = -1.0;
	this->npairs = 0;
	this->nex = -1;
	this->nsfg = -1;
	
	this->lwall = -1.0e6;
	this->lwall_cv = -1;
	this->lwall_n = 4;
	this->lwall_eps = 0.01;
	this->lwall_k = 5;

	this->use_restraint = 0;
	this->restraint_value = 0.0;
	this->restraint_n = 5;
	this->restraint_eps = 0.1;
	this->restraint_k = 40;
	this->restraint_pref = 0.0;
	
	this->useNN = -1;
	this->choose_exvar = -1;
	this->useMetadyn = 0;
	this->bias_stride = -1;
	this->bias_read = 0;
	this->bias_readfile = ' ';
	
	this->Rskappa = -1;
	this->eta = -1;
	this->nmol = -1;
	this->center = -1;
	this->natm = -1;
	this->nsfg2CO = -1;
	this->nsfg2NN = -1;
	this->nsfg3CO = -1;
	this->nsfg3NN = -1;
	this->nsfg2point = -1;
	this->nsfg3point = -1;
	this->triclinic = 0;
}

//Destructor
DAFED::CParameter::~CParameter()
{
}

//----------------------------Class CGgmt-----------------------------------------

//Constructor
DAFED::CGgmt::CGgmt()
{
	this->Q1 = 1.0;
	this->Q2 = 1.0;

	this->n_respa_ggmt = 1;
}

//Destructor
DAFED::CGgmt::~CGgmt()
{
}


//----------------------------Class CHisto----------------------------------------

//Constructor
DAFED::CHisto::CHisto()
{
	this->nbin = -1;
	this->nlat = -1;
	this->min = 0.0;
	this->max = 0.0;
	this->binwidth = -1.0;
	this->NconvDim = -1;
}

//Destructor
DAFED::CHisto::~CHisto()
{
}


//JR start: included for ramping s up/down
//--------------------------Class CRamp-------------------------------------------

//Constructor
DAFED::CRamp::CRamp()
{
	this->smin = 0.0;
	this->smax = 1.0;
	this->ds = 1e-3;
	this->steps = 10;
	this->steps0 = 10;
	this->count = 0;
	this->switchme = 1.0;
}

//Destructor
DAFED::CRamp::~CRamp()
{
}

//JR end: included for ramping s up/down


//----------------------------Class EXvariables-----------------------------------

//Constructor
DAFED::EXvariables::EXvariables()
{
	//this->neighbourdistance = -1.0;
	this->x = 0.0;
	this->v = 0.0;
	this->f = 0.0;
	this->fbias = 0.0;
	this->fcoup = 0.0;
	this->fwall = 0.0;
	this->fLwall = 0.0;
	this->fUwall = 0.0;
	this->fatom = 0.0;
	this->frestrain = 0.0;

	this->kappa = 0.0;
	this->gamma = 0.0;
	this->tau = 0.0;
	this->temp = 0.0;
	this->m = 0.0;
	this->c1_lgvn = 0.0;
	this->c2_lgvn = 0.0;

	this->n_respa = 1;

	this->lwall = -1.0e6;
	this->lwall_n = 4;
	this->lwall_eps = 0.01;
	this->lwall_k = 5;
	this->lwall_pref = 0.0;
	this->lwall_cv = -1;
	
	this->uwall = 1.0e6;
	this->uwall_n = 4;
	this->uwall_eps = 0.01;
	this->uwall_k = 5;
	this->uwall_pref = 0.0;
	this->uwall_cv = -1;

}

//Destructor
DAFED::EXvariables::~EXvariables()
{
}

//-------------------------CLASS QLMderiv------------------------------------------
//Constructor
//DAFED::QLMderiv::QLMderiv()
//Destructor

DAFED::NNvariables::NNvariables()
{}

DAFED::NNvariables::~NNvariables()
{}

DAFED::CNeighvariables::CNeighvariables()
{
}


DAFED::CNeighvariables::~CNeighvariables()
{
}

