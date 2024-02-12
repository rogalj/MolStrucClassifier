#include "readinput.h"


//--------------------------------------------------------------------------
// get number of atoms from lammps file
void get_natoms_lammps(const string &filename, int &natoms )
{
	int nop;
	char dummy_char[256];                //dummy line
	ifstream confFile;
	
	confFile.open(filename.c_str(),ifstream::in);

	if (confFile.is_open()){ 
		//the first 3 lines are dummy lines
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//the 4th line contains the number of particles
		confFile >> nop;
		natoms = nop;
	}
	else{
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
	confFile.close();
}








//--------------------------------------------------------------------------
// reading a lammps dump file as written in the TPS_wrapper
void readParticleFile_lammps(const string &filename, double **pos, int &natoms, double *box )
{
	double posx,posy,posz;
	//int nop;
	double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
	double dummy;                        //dummy variable
	char dummy_char[256];                //dummy line
	ifstream confFile;
	//int i,j;
	
	confFile.open(filename.c_str(),ifstream::in);


	if (confFile.is_open()){ 
		//the first 3 lines are dummy lines
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//the 4th line contains the number of particles
		confFile >> natoms;
		//initialize class for particles, etc.
		//this->initializeMolecules(nop);

		//again more dummy lines (also one extra after reading a value...
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//get the box, upper and lower value 
		confFile >> xsizeinf;
		confFile >> xsizesup;
		confFile >> ysizeinf;
		confFile >> ysizesup;
		confFile >> zsizeinf;
		confFile >> zsizesup;
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//cout<<xsizeinf<<"\n"<<xsizesup<<"\n"<<ysizeinf<<"\n"<<ysizesup<<"\n"<<zsizeinf<<"\n"<<zsizesup<<""<<"\n";  

		//determine box length
		box[0] = xsizesup - xsizeinf;
		box[1] = ysizesup - ysizeinf;
		box[2] = zsizesup - zsizeinf;

		//so lets read the particles positions
		for (int ti = 0;ti<natoms;ti++){
			//in the current format the first three values are the id, type and mass
			confFile>>dummy;
			confFile>>dummy;
			confFile>>dummy;
			//then comes x,y,z
			confFile>>posx;
			confFile>>posy;
			confFile>>posz;
			//then 3 more values for the velocities
			confFile>>dummy;
			confFile>>dummy;
			confFile>>dummy;
			//then 3 more for the forces (if printed)
			confFile>>dummy;
			confFile>>dummy;
			confFile>>dummy;


			//cout<<posx<<"\n"<<posy<<"\n"<<posz<<"\n";
     
			pos[ti][0] = posx;
			pos[ti][1] = posy;
			pos[ti][2] = posz;
		}
	}
	else{
		//cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}

//--------------------------------------------------------------------------
// reading a lammps dump file as written in the TPS_wrapper
void readParticleFile_lammps_aq6(const string &filename, double **pos, int &natoms, double *box, double *aq6 )
{
	double posx,posy,posz;
	//int nop;
	double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
	double dummy;                        //dummy variable
	char dummy_char[256];                //dummy line
	ifstream confFile;
	//int i,j;
	
	confFile.open(filename.c_str(),ifstream::in);


	if (confFile.is_open()){ 
		//the first 3 lines are dummy lines
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//the 4th line contains the number of particles
		confFile >> natoms;
		//initialize class for particles, etc.
		//this->initializeMolecules(nop);

		//again more dummy lines (also one extra after reading a value...
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//get the box, upper and lower value 
		confFile >> xsizeinf;
		confFile >> xsizesup;
		confFile >> ysizeinf;
		confFile >> ysizesup;
		confFile >> zsizeinf;
		confFile >> zsizesup;
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//cout<<xsizeinf<<"\n"<<xsizesup<<"\n"<<ysizeinf<<"\n"<<ysizesup<<"\n"<<zsizeinf<<"\n"<<zsizesup<<""<<"\n";  

		//determine box length
		box[0] = xsizesup - xsizeinf;
		box[1] = ysizesup - ysizeinf;
		box[2] = zsizesup - zsizeinf;

		//so lets read the particles positions
		for (int ti = 0;ti<natoms;ti++){
			//in the current format the first three values are the id, type and mass
			confFile>>dummy;
			confFile>>dummy;
			confFile>>dummy;
			//then comes x,y,z
			confFile>>posx;
			confFile>>posy;
			confFile>>posz;
			//then aq6
			confFile>>aq6[ti];


			//cout<<posx<<"\n"<<posy<<"\n"<<posz<<"\n";
     
			pos[ti][0] = posx;
			pos[ti][1] = posy;
			pos[ti][2] = posz;
		}
	}
	else{
		//cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}


void readParticleFile_sa_vec(const string &filename, double **pos, int &natoms, int *moltype, int *eletype, double *box, int &timeframe)
{
	double posx,posy,posz;
	int tempmol, tempele;
	//int nop;
	double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
	double dummy;                        //dummy variable
	char dummy_char[256];                //dummy line
        ifstream confFile;
        //int i,j;
 
	confFile.open(filename.c_str(),ifstream::in);
	
        if (confFile.is_open()){
  	      //the first 3 lines are dummy lines
              confFile.getline(dummy_char,256);
              confFile >> timeframe;
	      //confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              //the 4th line contains the number of particles
              confFile >> natoms;
              //initialize class for particles, etc.
              //this->initializeMolecules(nop);

              //again more dummy lines (also one extra after reading a value...
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              //get the box, upper and lower value 
              confFile >> xsizeinf;
              confFile >> xsizesup;
              confFile >> ysizeinf;
              confFile >> ysizesup;
              confFile >> zsizeinf;
              confFile >> zsizesup;
	      //cout << xsizesup;
              confFile.getline(dummy_char,256);
	      confFile.getline(dummy_char,256);
	      //cout<<xsizeinf<<"\n"<<xsizesup<<"\n"<<ysizeinf<<"\n"<<ysizesup<<"\n"<<zsizeinf<<"\n"<<zsizesup<<""<<"\n";  
	      
	      //determine box length
	      box[0] = xsizesup - xsizeinf;
	      box[1] = ysizesup - ysizeinf;
	      box[2] = zsizesup - zsizeinf;
	      //cout<<box[0]<<"\n"<<box[1]<<"\n"<<box[2]<<"\n";
	      //so lets read the particles positions
	      for (int ti = 0;ti<natoms;ti++){
	      //in the current format the first three values are the id, type and mass
	      	confFile>>dummy;
		confFile>>tempmol;
		confFile>>tempele;
	        //confFile>>dummy;
	      	//confFile>>dummy;
	        //then comes x,y,z
	        confFile>>posx;
	        confFile>>posy;
	        confFile>>posz;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;

	        //then aq6
	        //confFile>>vec[ti];
	        moltype[ti] = tempmol -1;
		eletype[ti] = tempele;
		//cout << moltype[ti] << "\n\n";	
		//cout << eletype[ti] << "\n";
		//cout<<posx<<"\n"<<posy<<"\n"<<posz<<"\n";
		pos[ti][0] = posx;
		pos[ti][1] = posy;
		pos[ti][2] = posz;
		}
	}
	else{
		//cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}



void readParticleFile_lammps_vec(const string &filename, double **pos, int &natoms, int *moltype, int *eletype, double *box)
{
	double posx,posy,posz;
	int tempmol, tempele;
	//int nop;
	double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
	double dummy;                        //dummy variable
	char dummy_char[256];                //dummy line
        ifstream confFile;
        //int i,j;
 
	confFile.open(filename.c_str(),ifstream::in);
	
        if (confFile.is_open()){
  	      //the first 3 lines are dummy lines
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              //the 4th line contains the number of particles
              confFile >> natoms;
              //initialize class for particles, etc.
              //this->initializeMolecules(nop);

              //again more dummy lines (also one extra after reading a value...
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              //get the box, upper and lower value 
              confFile >> xsizeinf;
              confFile >> xsizesup;
              confFile >> ysizeinf;
              confFile >> ysizesup;
              confFile >> zsizeinf;
              confFile >> zsizesup;
	      //cout << xsizesup;
              confFile.getline(dummy_char,256);
	      confFile.getline(dummy_char,256);
	      //cout<<xsizeinf<<"\n"<<xsizesup<<"\n"<<ysizeinf<<"\n"<<ysizesup<<"\n"<<zsizeinf<<"\n"<<zsizesup<<""<<"\n";  
	      
	      //determine box length
	      box[0] = xsizesup - xsizeinf;
	      box[1] = ysizesup - ysizeinf;
	      box[2] = zsizesup - zsizeinf;
	      //cout<<box[0]<<"\n"<<box[1]<<"\n"<<box[2]<<"\n";
	      //so lets read the particles positions
	      for (int ti = 0;ti<natoms;ti++){
	      //in the current format the first three values are the id, type and mass
	      	confFile>>dummy;
		confFile>>tempmol;
		confFile>>tempele;
	        //confFile>>dummy;
	      	//confFile>>dummy;
	        //then comes x,y,z
	        confFile>>posx;
	        confFile>>posy;
	        confFile>>posz;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;

	        //then aq6
	        //confFile>>vec[ti];
	        moltype[ti] = tempmol -1;
		eletype[ti] = tempele;
		//cout << moltype[ti] << "\n\n";	
		//cout << eletype[ti] << "\n";
		//cout<<posx<<"\n"<<posy<<"\n"<<posz<<"\n";
		pos[ti][0] = posx;
		pos[ti][1] = posy;
		pos[ti][2] = posz;
		}
	}
	else{
		//cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}


void readParticleFile_lammps_vec_tri(const string &filename, double **pos, int &natoms, int *moltype, int *eletype, double *box)
{
	double posx,posy,posz;
	int tempmol, tempele;
	//int nop;
	double xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,yz,xz,xlo,xhi,ylo,yhi,zlo,zhi;
	double dummy;                        //dummy variable
	char dummy_char[256];                //dummy line
        ifstream confFile;
 
	confFile.open(filename.c_str(),ifstream::in);
	
        if (confFile.is_open()){
  	      //the first 3 lines are dummy lines
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              //the 4th line contains the number of particles
              confFile >> natoms;
              //initialize class for particles, etc.
              //this->initializeMolecules(nop);

              //again more dummy lines (also one extra after reading a value...
              confFile.getline(dummy_char,256);
              confFile.getline(dummy_char,256);
              //get the box, upper and lower value 
              confFile >> xlo_bound;
              confFile >> xhi_bound;
              confFile >> xy;
              confFile >> ylo_bound;
              confFile >> yhi_bound;
              confFile >> xz;
              confFile >> zlo_bound;
              confFile >> zhi_bound;
              confFile >> yz;
	      xlo = xlo_bound-minx(xy,xz);
	      xhi = xhi_bound-maxx(xy,xz);
	      ylo = ylo_bound-miny(yz);
	      yhi = yhi_bound-maxy(yz);
	      zlo = zlo_bound;
	      zhi = zhi_bound;
              confFile.getline(dummy_char,256);
	      confFile.getline(dummy_char,256);
	      
	      //determine box length
	      box[0] = xhi - xlo;
	      box[1] = yhi - ylo;
	      box[2] = zhi - zlo;
	      box[3] = xy;
	      box[4] = xz;
	      box[5] = yz;
	      //so lets read the particles positions
	      for (int ti = 0;ti<natoms;ti++){
	      //in the current format the first three values are the id, type and mass
	      	confFile>>dummy;
		confFile>>tempmol;
		confFile>>tempele;
	        //then comes x,y,z
	        confFile>>posx;
	        confFile>>posy;
	        confFile>>posz;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;
		confFile>>dummy;

	        //then aq6
	        moltype[ti] = tempmol -1;
		eletype[ti] = tempele;
		pos[ti][0] = posx;
		pos[ti][1] = posy;
		pos[ti][2] = posz;
		}
	}
	else{
		//cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}


void readParticleFile_sa_vec_tri(const string &filename, double **pos, int &natoms, int *moltype, int *eletype, double *box, int& timeframe, int& timestep, double *boxinfo, int *G_ID)
{
	double posx,posy,posz;
	int tempid, tempmol, tempele;
	char dummy_char[256];                //dummy line
	double xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,yz,xz,xlo,xhi,ylo,yhi,zlo,zhi;
        ifstream confFile;
 
	confFile.open(filename.c_str(),ifstream::in);
	
	int i = 0;	
        if (confFile.is_open()){
	for (i = 0;i<timeframe*(natoms+9);i++)
	{
		confFile.getline(dummy_char,256);
	}
	confFile.getline(dummy_char,256);
	if(confFile.eof())
	{
		timeframe = -1;
		return;
	}
	confFile >> timestep;
	confFile.getline(dummy_char,256);
	confFile.getline(dummy_char,256);
	confFile >> natoms;
	confFile.getline(dummy_char,256);
	confFile.getline(dummy_char,256);
	confFile >> xlo_bound;
        confFile >> xhi_bound;
        confFile >> xy;
        confFile >> ylo_bound;
        confFile >> yhi_bound;
        confFile >> xz;
        confFile >> zlo_bound;
        confFile >> zhi_bound;
        confFile >> yz;
        xlo = xlo_bound-minx(xy,xz);
        xhi = xhi_bound-maxx(xy,xz);
        ylo = ylo_bound-miny(yz);
        yhi = yhi_bound-maxy(yz);
        zlo = zlo_bound;
        zhi = zhi_bound;
	box[0] = xhi - xlo;
	box[1] = yhi - ylo;
	box[2] = zhi - zlo;
	box[3] = xy;
	box[4] = xz;
	box[5] = yz;
	       	
	boxinfo[0] = xlo_bound;
	boxinfo[1] = xhi_bound;
	boxinfo[2] = xy;
	boxinfo[3] = ylo_bound;
	boxinfo[4] = yhi_bound;
	boxinfo[5] = xz;
	boxinfo[6] = zlo_bound;
	boxinfo[7] = zhi_bound;
	boxinfo[8] = yz;
	confFile.getline(dummy_char,256);
	confFile.getline(dummy_char,256);
	
	      for (int ti = 0;ti<natoms;ti++){
	      //in the current format the first three values are the id, type and mass
	      	confFile>>tempid;               //Global ID;
		confFile>>tempmol;
		confFile>>tempele;
	        //then comes x,y,z
	        confFile>>posx;
	        confFile>>posy;
	        confFile>>posz;
		confFile.getline(dummy_char,256);
		
		G_ID[ti] = tempid-1;	
	        moltype[ti] = tempmol -1;
		eletype[ti] = tempele;
		pos[ti][0] = posx;
		pos[ti][1] = posy;
		pos[ti][2] = posz;
		}
	
	timeframe++;
	}
	else{
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}



void readParticleFile_sa_vec_ortho(const string &filename, double **pos, int &natoms, int *moltype, int *eletype, double *box, int& timeframe, int& timestep, double *boxinfo, int *G_ID)
{
	double posx,posy,posz;
	int tempid, tempmol, tempele;
	char dummy_char[256];                //dummy line
        double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
        ifstream confFile;
 
	confFile.open(filename.c_str(),ifstream::in);
	
	int i = 0;	
        if (confFile.is_open()){
	for (i = 0;i<timeframe*(natoms+9);i++)
	{
		confFile.getline(dummy_char,256);
	}
	confFile.getline(dummy_char,256);
	if(confFile.eof())
	{
		timeframe = -1;
		return;
	}
	confFile >> timestep;
	confFile.getline(dummy_char,256);
	confFile.getline(dummy_char,256);
	confFile >> natoms;
	confFile.getline(dummy_char,256);
	confFile.getline(dummy_char,256);
              confFile >> xsizeinf;
              confFile >> xsizesup;
              confFile >> ysizeinf;
              confFile >> ysizesup;
              confFile >> zsizeinf;
              confFile >> zsizesup;
	      box[0] = xsizesup - xsizeinf;
	      box[1] = ysizesup - ysizeinf;
	      box[2] = zsizesup - zsizeinf;
        
	boxinfo[0] = xsizeinf;
	boxinfo[1] = xsizesup;
	boxinfo[2] = ysizeinf;
	boxinfo[3] = ysizesup;
	boxinfo[4] = zsizeinf;
	boxinfo[5] = zsizesup;
	confFile.getline(dummy_char,256);
	confFile.getline(dummy_char,256);
	
	      for (int ti = 0;ti<natoms;ti++){
	      //in the current format the first three values are the id, type and mass
	      	confFile>>tempid;		//Global ID
		confFile>>tempmol;
		confFile>>tempele;
	        //then comes x,y,z
	        confFile>>posx;
	        confFile>>posy;
	        confFile>>posz;
		
		confFile.getline(dummy_char,256);
		G_ID[ti] = tempid-1;
	        moltype[ti] = tempmol -1;
		eletype[ti] = tempele;
		pos[ti][0] = posx;
		pos[ti][1] = posy;
		pos[ti][2] = posz;
		}
	
	timeframe++;
	}
	else{
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}


//--------------------------------------------------------------------------
// reading a lammps dump file with 22 q values
void readParticleFile_lammps_allq(const string &filename, double **pos, int &natoms, double *box, double **qval )
{
	double posx,posy,posz;
	//int nop;
	double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
	double dummy;                        //dummy variable
	char dummy_char[256];                //dummy line
	ifstream confFile;
	//int i,j;
	
	confFile.open(filename.c_str(),ifstream::in);


	if (confFile.is_open()){ 
		//the first 3 lines are dummy lines
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//the 4th line contains the number of particles
		confFile >> natoms;
		//initialize class for particles, etc.
		//this->initializeMolecules(nop);

		//again more dummy lines (also one extra after reading a value...
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//get the box, upper and lower value 
		confFile >> xsizeinf;
		confFile >> xsizesup;
		confFile >> ysizeinf;
		confFile >> ysizesup;
		confFile >> zsizeinf;
		confFile >> zsizesup;
		confFile.getline(dummy_char,256);
		confFile.getline(dummy_char,256);
		//cout<<xsizeinf<<"\n"<<xsizesup<<"\n"<<ysizeinf<<"\n"<<ysizesup<<"\n"<<zsizeinf<<"\n"<<zsizesup<<""<<"\n";  

		//determine box length
		box[0] = xsizesup - xsizeinf;
		box[1] = ysizesup - ysizeinf;
		box[2] = zsizesup - zsizeinf;

		//so lets read the particles positions
		for (int ti = 0;ti<natoms;ti++){
			//in the current format the first three values are the id, type and mass
			confFile>>dummy;
			confFile>>dummy;
			confFile>>dummy;
			//then comes x,y,z
			confFile>>posx;
			confFile>>posy;
			confFile>>posz;
			//then 22 q values
			for(int j=0 ; j<22; j++){
				confFile>>qval[ti][j];
			}


			//cout<<posx<<"\n"<<posy<<"\n"<<posz<<"\n";
     
			pos[ti][0] = posx;
			pos[ti][1] = posy;
			pos[ti][2] = posz;
		}
	}
	else{
		//cerr << "Fatal Error : cannot open the file " <<  this->parameter->xyzFile << "\n";
		cerr << "Fatal Error : cannot open the file " <<  filename << "\n";
		cerr << "Exiting programme..."<<endl;
		exit(1);
	}
}

double minx(double xy, double xz)
{
	if(xy <0  && xz <0)
		return xy + xz;
	else if(xy <0|| xz <0)
	{
		if (xy < xz)
			return xy;
		else 
			return xz;
	}
	else 
		return 0;
}
double maxx(double xy, double xz)
{
	if(xy >0  && xz >0)
		return xy + xz;
	else if(xy >0|| xz >0)
	{
		if (xy > xz)
			return xy;
		else 
			return xz;
	}
	else 
		return 0;
}

double miny(double yz)
{
	if (yz <0)
		return yz;
	else 
		return 0;
}
double maxy(double yz)
{
	if (yz >0)
		return yz;
	else 
		return 0;
}
