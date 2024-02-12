//Need both header and dafed. 
#include "header.h"	
#include "dafed.h"


//This is mimicing lammps dot3 and len3. Just making a dot product of two vector(3D) and length of the vector
double dot3(double v1[],double v2[]){
	        return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}
double len3(double v[]){
	        return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

//from here up to "//------------------------------------------------------------------------------------------------------------------------------------------------------------------" are proto code
//This function will calculate G2 for CO. Data will be stored in molecules as molec[ti].G2vecCO
void DAFED::symmetryfunc_g2_fc0_vecCO_lmp(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter, int *COvectype, CNeighvariables &lmpneigh){
        int mi,mj, jmol,tj;			//mi for loop all the molec, mj for loop over all neighbor molecule and jmol for directing right neighbor molecule 
        double cosRs = parameter.Rskappa;	//cos(Rs)
        double eta = parameter.eta;		//eta 
        double g2;				//sum of G2 for the target molecule
        double dr;				//dr is cosvec1-cosRs
	double new_molvec1[3];     		//vector of the target CO
        double new_molvec2[3];     		//vector of the neighbor CO
        double cosvec1;				//cos(theta_ikjl)
        int ii;				//center of the targe and neighbor molecule 
        for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
                ii = molec[mi].center;		//center of the target molecule.
                molec[mi].G2vecCO = 0;		//No need to have, but this is used in here to make sure G2becCO is 0.
                if (molec[mi].center==-1)break;
		for (int k = 0; k < 3; k++)	//xyz component of the vector. 
                {
	                new_molvec1[k] = molecules[molec[mi].center].COvec[k];		//assigning the molvec
                }
                g2 = 0;				//g2 start from 0 and adds up in the g2
                for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			tj =jmol;
			if (ii==tj)				//if it is the same molecule, skip. In principle, this would be never called.
				continue;
			for (int k = 0; k < 3; k++)		//xyz component of the vector for neighbor
			{
	                	new_molvec2[k] = molecules[tj].COvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));		//cos(theta_ikjl)
			dr = cosvec1-cosRs;									//dr is cosvec1-cosRs
			g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];						//sum of G2 for the target molecule
		}
		molec[mi].G2vecCO = g2;										//assigning sum of the g2 to each molecule 
	}
}

void DAFED::symmetryfunc_g2_fc0_vecNN_lmp(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter, int *NNvectype,CNeighvariables &lmpneigh){
//add type in the molecule
//add id in molecule
//think about what to do with ichunk and cvecchunk. foe cvecchunk, think about pu_molvecall. Maybe get type from other places?
//add mask on molecule?

	int mi,mj, jmol, tj;			//mi for loop all the molec, mj for loop over all neighbor molecule and jmol for directing right neighbor molecule
	double cosRs = parameter.Rskappa;	//cos(Rs)
	double eta = parameter.eta;		//eta from parameter
	double g2;				//sum of G2 for the target molecule
	double dr;				//dr is cosvec1-cosRs
        double new_molvec1[3];     		//vector of the target NN
        double new_molvec2[3];     		//vector of the neighbor NN
        double cosvec1;				//cos(theta_ikjl)
	int ii;				//center of the targe and neighbor molecule
	for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
		ii = molec[mi].center;		//center of the target molecule.
                if (molec[mi].center==-1)break;
		for (int k = 0; k < 3; k++)	//xyz component of the vector.
		{
	                new_molvec1[k] = molecules[molec[mi].center].NNvec[k];		//assigning the molvec
		}
		g2 = 0;						//g2 start from 0 and adds up in the g2
		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules. 
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			tj = jmol;
			if (ii==tj)				//if it is the same molecule, skip. In principle, this would be never called.
				continue;
			for (int k = 0; k < 3; k++)		//xyz component of the vector for neighbor
			{
	                	new_molvec2[k] = molecules[tj].NNvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));		//cos(theta_ikjl)
			dr = cosvec1-cosRs;									//dr is cosvec1-cosRs
			g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];						//sum of G2 for the target molecule
		}
		molec[mi].G2vecNN = g2;										//assigning sum of the g2 to each molecule
	}
}


//This is for point representation. Very similar to one above or one Jutta wrote. 
void DAFED::symfg2_lmp(CAtom *atoms, Cmolpoint *molec, int &natoms, CParameter &parameter,CNeighvariables &lmpneigh){
	int mi,mj, jmol,tj;
	double Rs = parameter.Rskappa;
	double eta = parameter.eta;
	double g2;
	double dr;
	int ii;
	for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
		ii = molec[mi].center;
		g2 = 0;
                if (molec[mi].center==-1)break;
		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.     		
			jmol = molec[mi].neighbors[mj];
			
			tj = jmol;
			if (ii==tj)
				continue;
			dr = molec[mi].neighdist[mj]-Rs;
			g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];
		}
		molec[mi].point2 = g2;
	}
}
void DAFED::symfg3_lmp(CAtom *atoms, Cmolpoint *molec, int &natoms, CParameter &parameter,CNeighvariables &lmpneigh){
	int mi,mj, jmol,tj;
	double kappa = parameter.Rskappa;
	double g3;
	double r;
	int ii;
	for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
		ii = molec[mi].center;
		g3 = 0;
                if (molec[mi].center==-1)break;
		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];		
			tj = jmol;
			if (ii==tj)
				continue;
			
			r = abs(molec[mi].neighdist[mj]);
			g3 += cos(kappa*r)*molec[mi].fcv[mj];
		}
		molec[mi].point3 = g3;		
	}
}


void DAFED::symmetryfunc_g3_fc0_vecNN_lmp(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter, int *NNvectype,CNeighvariables &lmpneigh){
	int mi,mj, jmol,tj;
	double kappa = parameter.Rskappa;
	double g3;
        double new_molvec1[3];     // = cvecchunk1->pu_molvecall;     make it somewhere in the loop
        double new_molvec2[3];     // = cvecchunk2->pu_molvecall;
        double cosvec1;
	int ii;
	for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
		ii = molec[mi].center;
                if (molec[mi].center==-1)break;
		for (int k = 0; k < 3; k++)
		{
	                new_molvec1[k] = molecules[molec[mi].center].NNvec[k];		//assigning the molvec
		}
		g3 = 0;
		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];
			tj = jmol;
			if (ii==tj)
				continue;
			for (int k = 0; k < 3; k++)
			{
	                	new_molvec2[k] = molecules[tj].NNvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));
			g3 += cos(kappa*cosvec1)*molec[mi].fcv[mj];
		}
		molec[mi].G3vecNN = g3;
	}
}


void DAFED::symmetryfunc_g3_fc0_vecCO_lmp(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter, int *COvectype,CNeighvariables &lmpneigh){
        int ti,tj;
        double kappa = parameter.Rskappa;
        double g3;
        double new_molvec1[3];     // = cvecchunk1->pu_molvecall;     make it somewhere in the loop
        double new_molvec2[3];     // = cvecchunk2->pu_molvecall;
        double cosvec1;
        int ii,jmol;
        for(ti=0;ti<lmpneigh.lmpmol;ti++){       //loop over all molecules
                ii = molec[ti].center;
                if (molec[ti].center==-1)break;
		molec[ti].G3vecCO = 0;
		for (int k = 0; k < 3; k++)
                {
	                new_molvec1[k] = molecules[molec[ti].center].COvec[k];		//assigning the molvec
                }
                g3 = 0;
                for (tj = 0; tj < molec[ti].n_neighbors;tj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[ti].neighbors[tj];
			if (ii==jmol)
				continue;
			for (int k = 0; k < 3; k++)
			{
	                	new_molvec2[k] = molecules[jmol].COvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));
			g3 += cos(kappa*cosvec1)*molec[ti].fcv[tj];
		}
		molec[ti].G3vecCO = g3;
	}
}


void DAFED::get_total_derivatives_molvec_lmp(CAtom *atoms, Cmolpoint *molec, int &natoms,CParameter &parameter, int *vectypeCO, int *vectypeNN, EXvariables *exvar, double *RskappaLst, double *etaLst, CNeighvariables &lmpneigh){
	int mi,mj,mk;
	int jmol;
        double Rskappa;
        double eta;
	double dsfg;
	double dr;
	double new_molvec1CO[3];
	double new_molvec1NN[3];
	double new_molvec2CO[3];
	double new_molvec2NN[3];
	double CCO,CNN, dC, ACO,ANN, dA;
        int atmi = 0,atmj = 0,atmk = 0,atmm = 0,atmo =0;
        double rikjl,rmonp, magrik, magrjl, magrmo, magrnp,nume1CO,nume1NN,nume2CO,nume2NN, denoCO,denoNN;
	int nsfg = parameter.nsfg;
	int nsfg2point = parameter.nsfg2point;
	int nsfg2CO = parameter.nsfg2CO;
	int nsfg2NN = parameter.nsfg2NN;
	int nsfg3point = parameter.nsfg3point;
	int nsfg3CO = parameter.nsfg3CO;
	int nsfg3NN = parameter.nsfg3NN;
	double qII,qJJ;
	int outi;
	int outNN_index[2];
	double dgdxI[4][3];
	double dgdxJ[4][3];
	double nmol = parameter.nmol;
	double lmpmol = lmpneigh.lmpmol;
	outNN_index[0] =0;		//output of urea I in the NN
	outNN_index[1] =1;		//output of urea IV in the NN
	for(int iex=0;iex<parameter.nex;iex++){
	        for(int ti=0;ti<lmpneigh.nall;ti++){
			for(int i=0;i<3;i++){
				exvar[iex].dQ[ti][i] = 0.0;
			}
		}
	}
	for(mi=0;mi<lmpmol;mi++)		//loop over all molecules original was parameter.nmol, but changed to 10 for testing.
	{
		if (molec[mi].center ==-1)break;
		
		atmi = molec[mi].atom_ID[0];
		atmk = molec[mi].atom_ID[1];
		atmm = molec[mi].atom_ID[2];
		atmo = molec[mi].atom_ID[3];
		for (int k = 0; k < 3; k++)		//get molvector of target molecule
		{             
	                new_molvec1CO[k] = atoms[atmi].COvec[k];		//assigning the molvec
	                new_molvec1NN[k] = atoms[atmi].NNvec[k];		//assigning the molvec
		}
		for (mj = 0; mj < molec[mi].n_neighbors;mj++)		//loop over all the neighbor molecules.
		{ 
	                jmol = molec[mi].neighbors[mj];////////////////////////////	//actually, this is not mol. this is atmj(carbon atom only)
			atmj = jmol;
			if (atmi==atmj)                     //Skip if target and neighbor is the same molecule
				continue;	
			for (int k = 0; k < 3; k++)
			{
	                	new_molvec2CO[k] = atoms[atmj].COvec[k];		//assigning the molvec
	                	new_molvec2NN[k] = atoms[atmj].NNvec[k];		//assigning the molvec
			}
			rikjl = dot3(new_molvec1CO,new_molvec2CO);
			rmonp = dot3(new_molvec1NN,new_molvec2NN);
			magrik = (len3(new_molvec1CO));
			magrjl = (len3(new_molvec2CO));
			magrmo = (len3(new_molvec1NN));
			magrnp = (len3(new_molvec2NN));

			nume1CO = magrik*magrjl;
			nume2CO = rikjl*magrjl/magrik;
			
			nume1NN = magrmo*magrnp;
			nume2NN = rmonp*magrnp/magrmo;
			
			denoCO = magrik*magrjl*magrik*magrjl;
			denoNN = magrmo*magrnp*magrmo*magrnp;
			for (mk = 0; mk < nsfg; mk++)
			{
			
				for (int ti = 0; ti < 4; ti++)
				{
					for(int k = 0;k<3;k++)  //for x,y,z
					{
						dgdxI[ti][k] = 0;
						dgdxJ[ti][k] = 0;	
					}
				}
				Rskappa = RskappaLst[mk];
	                        eta = etaLst[mk];
				if (mk < nsfg2point)		//point g2
				{
					dr = molec[mi].neighdist[mj]-Rskappa;
					for (int k = 0; k < 3; k++)
					{
						dsfg = exp(-eta*dr*dr)*(2.0*eta*dr*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO)	//CO g2
				{

					for (int k = 0; k < 3; k++)
					{
						ACO = rikjl/(magrik*magrjl) - Rskappa;
						dA = (-nume1CO*new_molvec2CO[k] + nume2CO*new_molvec1CO[k])/denoCO;					//derivative of A 
						dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);		//derivative of g by xi,yi,zi
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
						dA = -dA;								//When we have a derivative of G2I or G2J with atmk, dA will be just flipped. 
                                		dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);					//Same calc but dA sign is opposite
						dgdxI[1][k] = dsfg;
						dgdxJ[1][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN)	//NN g2
				{
				
					for (int k = 0; k < 3; k++)
					{
						ANN = rmonp/(magrmo*magrnp) - Rskappa;
						dA = (-nume1NN*new_molvec2NN[k] + nume2NN*new_molvec1NN[k])/denoNN;
        		                        dsfg = exp(-eta*ANN*ANN)*(molec[mi].dfcvdx[mj][k]);	////////check where i put negative sign.<-- I need to put negative on dfcvdx on dGJI part.
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxI[2][k] = dsfg;
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxJ[2][k] = dsfg;

						//For N2
						dA = -dA;
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxI[3][k] = dsfg;
        	                        	
						//dGJI/dxi for N2
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxJ[3][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN + nsfg3point)	//point g3
				{
					for (int k = 0; k < 3; k++)
					{
						dsfg = Rskappa*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj]*sin(Rskappa*molec[mi].neighdist[mj]) 
							+ cos(Rskappa*molec[mi].neighdist[mj])*molec[mi].dfcvdx[mj][k];
						dgdxI[0][k] = dsfg;
						dsfg =  dsfg;
						dgdxJ[0][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN + nsfg3point + nsfg3CO)	//CO g3
				{
						
					for (int k = 0; k < 3; k++)
					{
						CCO = Rskappa * rikjl/(magrik*magrjl);
						dC = Rskappa*(-nume1CO*new_molvec2CO[k] + nume2CO*new_molvec1CO[k])/denoCO;
			                        dsfg = -dC*sin(CCO)*molec[mi].fcv[mj] + cos(CCO)*molec[mi].dfcvdx[mj][k];
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
						dC = -dC;
						dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
						dgdxI[1][k] = dsfg;
						dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
						dgdxJ[1][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN + nsfg3point + nsfg3CO + nsfg3NN)	//NN g3
				{						
					for (int k = 0; k < 3; k++)
					{
						CNN = Rskappa * rmonp/(magrmo*magrnp);
						dC = Rskappa*(-nume1NN*new_molvec2NN[k] + nume2NN*new_molvec1NN[k])/denoNN;
	        	                        dsfg = cos(CNN)*molec[mi].dfcvdx[mj][k];
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
	                	                dsfg = -dC*sin(CNN)*molec[mi].fcv[mj];
						dgdxI[2][k] = dsfg;
						dgdxJ[2][k] = dsfg;
				//For N2
						dC = -dC;
						dsfg = -dC*sin(CNN)*molec[mi].fcv[mj];
						dgdxI[3][k] = dsfg;
						dsfg = -dC*sin(CNN)*molec[mi].fcv[mj];
						dgdxJ[3][k] = dsfg;
					}
				}
				else
					cout << "ERROR, exceeded the limit" <<endl;
				for (int iex = 0; iex < parameter.nex;iex++)
				{
					for (int k = 0; k < 3; k++)
					{
						outi = outNN_index[iex]*parameter.nsfg + mk;
						qII = atoms[atmi].NNgrad[outi]/nmol;
						qJJ = atoms[atmj].NNgrad[outi]/nmol;
						
						exvar[iex].dQ[atmi][k]+=qII*dgdxI[0][k]+qJJ*dgdxJ[0][k];
						exvar[iex].dQ[atmk][k]+=qII*dgdxI[1][k]+qJJ*dgdxJ[1][k];
						exvar[iex].dQ[atmm][k]+=qII*dgdxI[2][k]+qJJ*dgdxJ[2][k];
						exvar[iex].dQ[atmo][k]+=qII*dgdxI[3][k]+qJJ*dgdxJ[3][k];
					
						for(int ir=0;ir<3;ir++)
						{        //loop over x,y,z for position, add up virial
							exvar[iex].dvirial[ir][k] += (qII*(dgdxI[0][k]+dgdxI[1][k]+dgdxI[2][k]+dgdxI[3][k])
								+qJJ*(dgdxJ[0][k]+dgdxJ[1][k]+dgdxJ[2][k]+dgdxJ[3][k]))*molec[mi].diff[mj][ir]*0.5;
                                        	}//loop over xyz for virial
					}//loop over xyz for calc dQ
				}//loop over extended variables
			}//loop over all symmetry functions
		}//loop over neighbor molecules
	}//loop over target molecules
}









//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//stand alone code calculate each symmetry functions here up to "//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------" 
void DAFED::symmetryfunc_g2_fc0_vecCO_sa(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter){
        int mi,mj, jmol,tj;			//mi for loop all the molec, mj for loop over all neighbor molecule and jmol for directing right neighbor molecule 
        double cosRs = parameter.Rskappa;	//cos(Rs)
        double eta = parameter.eta;		//eta 
        double g2;				//sum of G2 for the target molecule
        double dr;				//dr is cosvec1-cosRs
	double new_molvec1[3];     		//vector of the target CO
        double new_molvec2[3];     		//vector of the neighbor CO
        double cosvec1;				//cos(theta_ikjl)
        int ii;				//center of the targe and neighbor molecule 
        for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
                ii = molec[mi].center;		//center of the target molecule.
                molec[mi].G2vecCO = 0;		//No need to have, but this is used in here to make sure G2becCO is 0.
                if (molec[mi].center==-1)break;
		for (int k = 0; k < 3; k++)	//xyz component of the vector. 
                {
	                new_molvec1[k] = molec[mi].COvec[k];		//assigning the molvec
                }
                g2 = 0;				//g2 start from 0 and adds up in the g2
                for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			tj =molec[jmol].center;			//center of the molec id.
			if (ii==tj)				//if it is the same molecule, skip. In principle, this would be never called.
				continue;
			for (int k = 0; k < 3; k++)		//xyz component of the vector for neighbor
			{
	                	new_molvec2[k] = molec[jmol].COvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));		//cos(theta_ikjl)
			dr = cosvec1-cosRs;									//dr is cosvec1-cosRs
			g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];						//sum of G2 for the target molecule
		}
		molec[mi].G2vecCO = g2;										//assigning sum of the g2 to each molecule 
	}
}

void DAFED::symmetryfunc_g2_fc0_vecNN_sa(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter){
        int mi,mj, jmol,tj;			//mi for loop all the molec, mj for loop over all neighbor molecule and jmol for directing right neighbor molecule 
        double cosRs = parameter.Rskappa;	//cos(Rs)
        double eta = parameter.eta;		//eta 
        double g2;				//sum of G2 for the target molecule
        double dr;				//dr is cosvec1-cosRs
	double new_molvec1[3];     		//vector of the target NN
        double new_molvec2[3];     		//vector of the neighbor NN
        double cosvec1;				//cos(theta_ikjl)
        int ii;				//center of the targe and neighbor molecule 
        for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
                ii = molec[mi].center;		//center of the target molecule.
                molec[mi].G2vecNN = 0;		//No need to have, but this is used in here to make sure G2vecNN is 0.
                if (molec[mi].center==-1)break;
		for (int k = 0; k < 3; k++)	//xyz component of the vector. 
                {
	                new_molvec1[k] = molec[mi].NNvec[k];		//assigning the molvec
                }
                g2 = 0;				//g2 start from 0 and adds up in the g2
                for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			tj =molec[jmol].center;			//center of the molec id.
			if (ii==tj)				//if it is the same molecule, skip. In principle, this would be never called.
				continue;
			for (int k = 0; k < 3; k++)		//xyz component of the vector for neighbor
			{
	                	new_molvec2[k] = molec[jmol].NNvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));		//cos(theta_ikjl)
			dr = cosvec1-cosRs;									//dr is cosvec1-cosRs
			g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];						//sum of G2 for the target molecule
		}
		molec[mi].G2vecNN = g2;										//assigning sum of the g2 to each molecule 
	}
}


void DAFED::symfg2_sa(CAtom *atoms, Cmolpoint *molec, int &natoms, CParameter &parameter){
	int mi,mj, jmol,tj;
	double Rs = parameter.Rskappa;
	double eta = parameter.eta;
	double g2;
	double dr;
	int ii;
	for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
		ii = molec[mi].center;
		g2 = 0;
                if (molec[mi].center==-1)break;
                for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			tj =molec[jmol].center;			//center of the molec id.
			
			if (ii==tj)
				continue;
			dr = molec[mi].neighdist[mj]-Rs;
			g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];
		}
		molec[mi].point2 = g2;
	}
}

void DAFED::symfg3_sa(CAtom *atoms, Cmolpoint *molec, int &natoms, CParameter &parameter){
	int mi,mj, jmol,tj;
	double kappa = parameter.Rskappa;
	double g3;
	double r;
	int ii;
	for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
		ii = molec[mi].center;
		g3 = 0;
                if (molec[mi].center==-1)break;
		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];		
			tj = molec[jmol].center;
			if (ii==tj)
				continue;
			
			r = abs(molec[mi].neighdist[mj]);
			g3 += cos(kappa*r)*molec[mi].fcv[mj];
		}
		molec[mi].point3 = g3;		
	}
}

void DAFED::symmetryfunc_g3_fc0_vecNN_sa(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter){
	int mi,mj, jmol,tj;
	double kappa = parameter.Rskappa;
	double g3;
        double new_molvec1[3];     // = cvecchunk1->pu_molvecall;     make it somewhere in the loop
        double new_molvec2[3];     // = cvecchunk2->pu_molvecall;
        double cosvec1;
	int ii;
	for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
		ii = molec[mi].center;
                if (molec[mi].center==-1)break;
		for (int k = 0; k < 3; k++)
		{
	                new_molvec1[k] = molec[mi].NNvec[k];		//assigning the molvec
		}
		g3 = 0;
		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];
			tj =molec[jmol].center;			//center of the molec id.
			if (ii==tj)
				continue;
			for (int k = 0; k < 3; k++)
			{
	                	new_molvec2[k] = molec[jmol].NNvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));
			g3 += cos(kappa*cosvec1)*molec[mi].fcv[mj];
		}
		molec[mi].G3vecNN = g3;
	}
}

void DAFED::symmetryfunc_g3_fc0_vecCO_sa(CAtom *molecules, Cmolpoint *molec, int &natoms, CParameter &parameter){
        int mi,mj,tj;
        double kappa = parameter.Rskappa;
        double g3;
        double new_molvec1[3];     // = cvecchunk1->pu_molvecall;     make it somewhere in the loop
        double new_molvec2[3];     // = cvecchunk2->pu_molvecall;
        double cosvec1;
        int ii,jmol;
        for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
                ii = molec[mi].center;
                if (molec[mi].center==-1)break;
		molec[mi].G3vecCO = 0;
		for (int k = 0; k < 3; k++)
                {
	                new_molvec1[k] = molec[mi].COvec[k];		//assigning the molvec
                }
                g3 = 0;
                for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
			jmol = molec[mi].neighbors[mj];
			tj =molec[jmol].center;			//center of the molec id.
			if (ii==tj)
				continue;
			for (int k = 0; k < 3; k++)
			{
	                	new_molvec2[k] = molec[jmol].COvec[k];		//assigning the molvec
			}
			cosvec1 = dot3(new_molvec1,new_molvec2)/(len3(new_molvec1)*len3(new_molvec2));
			g3 += cos(kappa*cosvec1)*molec[mi].fcv[mj];
		}
		molec[mi].G3vecCO = g3;
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//here we calculate the derivative of each symmetry function. 
void DAFED::get_total_derivatives_molvec_sa(CAtom *atoms, Cmolpoint *molec, int &natoms,CParameter &parameter, EXvariables *exvar, double *RskappaLst, double *etaLst){
	int mi,mj,mk;
	int jmol;
        double Rskappa;
        double eta;
	double dsfg;
	double dr;
	double new_molvec1CO[3];
	double new_molvec1NN[3];
	double new_molvec2CO[3];
	double new_molvec2NN[3];
	double CCO,CNN, dC, ACO,ANN, dA;
        int atmi = 0,atmk = 0,atmm = 0,atmo =0;
        double rikjl,rmonp, magrik, magrjl, magrmo, magrnp,nume1CO,nume1NN,nume2CO,nume2NN, denoCO,denoNN;
	int nsfg = parameter.nsfg;
	int nsfg2point = parameter.nsfg2point;
	int nsfg2CO = parameter.nsfg2CO;
	int nsfg2NN = parameter.nsfg2NN;
	int nsfg3point = parameter.nsfg3point;
	int nsfg3CO = parameter.nsfg3CO;
	int nsfg3NN = parameter.nsfg3NN;
	double qII,qJJ;
	int outi;
	int outNN_index[2];
	double dgdxI[4][3];
	double dgdxJ[4][3];
	double nmol = parameter.nmol;
	outNN_index[0] =0;		//output of urea I in the NN
	outNN_index[1] =1;		//output of urea IV in the NN
	for(int iex=0;iex<parameter.nex;iex++){
	        for(int ti=0;ti<natoms;ti++){	//	******This code does not work on parallel because of this line. If you want to run in parallel, change this to local + ghost ***********
			for(int i=0;i<3;i++){
				exvar[iex].dQ[ti][i] = 0.0;
			}
		}
	}
	for(mi=0;mi<nmol;mi++)		//loop over all molecules original was parameter.nmol, but changed to 10 for testing.
	{
		if (molec[mi].center ==-1)break;
		
		atmi = molec[mi].atom_ID[0];
		atmk = molec[mi].atom_ID[1];
		atmm = molec[mi].atom_ID[2];
		atmo = molec[mi].atom_ID[3];
		int icenter = molec[mi].center;
		for (int k = 0; k < 3; k++)		//get molvector of target molecule
		{             
	                new_molvec1CO[k] = atoms[icenter].COvec[k];		//assigning the molvec
	                new_molvec1NN[k] = atoms[icenter].NNvec[k];		//assigning the molvec
		}
		for (mj = 0; mj < molec[mi].n_neighbors;mj++)		//loop over all the neighbor molecules.
		{ 
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			int jcenter = molec[jmol].center;                 //center of the molec jmol id. 
			if (icenter==jcenter)                     //Skip if target and neighbor is the same molecule
				continue;	
			for (int k = 0; k < 3; k++)
			{
	                	new_molvec2CO[k] = atoms[jcenter].COvec[k];		//assigning the molvec
	                	new_molvec2NN[k] = atoms[jcenter].NNvec[k];		//assigning the molvec
			}
			rikjl = dot3(new_molvec1CO,new_molvec2CO);
			rmonp = dot3(new_molvec1NN,new_molvec2NN);
			magrik = (len3(new_molvec1CO));
			magrjl = (len3(new_molvec2CO));
			magrmo = (len3(new_molvec1NN));
			magrnp = (len3(new_molvec2NN));

			nume1CO = magrik*magrjl;
			nume2CO = rikjl*magrjl/magrik;
			
			nume1NN = magrmo*magrnp;
			nume2NN = rmonp*magrnp/magrmo;
			
			denoCO = magrik*magrjl*magrik*magrjl;
			denoNN = magrmo*magrnp*magrmo*magrnp;
			for (mk = 0; mk < nsfg; mk++)
			{
			
				for (int ti = 0; ti < 4; ti++)
				{
					for(int k = 0;k<3;k++)  //for x,y,z
					{
						dgdxI[ti][k] = 0;
						dgdxJ[ti][k] = 0;	
					}
				}
				Rskappa = RskappaLst[mk];
	                        eta = etaLst[mk];
				if (mk < nsfg2point)		//point g2
				{
					dr = molec[mi].neighdist[mj]-Rskappa;
					for (int k = 0; k < 3; k++)
					{
						dsfg = exp(-eta*dr*dr)*(2.0*eta*dr*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO)	//CO g2
				{

					for (int k = 0; k < 3; k++)
					{
						ACO = rikjl/(magrik*magrjl) - Rskappa;
						dA = (-nume1CO*new_molvec2CO[k] + nume2CO*new_molvec1CO[k])/denoCO;					//derivative of A 
						dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);		//derivative of g by xi,yi,zi
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
						dA = -dA;								//When we have a derivative of G2I or G2J with atmk, dA will be just flipped. 
                                		dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);					//Same calc but dA sign is opposite
						dgdxI[1][k] = dsfg;
						dgdxJ[1][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN)	//NN g2
				{
				
					for (int k = 0; k < 3; k++)
					{
						ANN = rmonp/(magrmo*magrnp) - Rskappa;
						dA = (-nume1NN*new_molvec2NN[k] + nume2NN*new_molvec1NN[k])/denoNN;
        		                        dsfg = exp(-eta*ANN*ANN)*(molec[mi].dfcvdx[mj][k]);	////////check where i put negative sign.<-- I need to put negative on dfcvdx on dGJI part.
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxI[2][k] = dsfg;
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxJ[2][k] = dsfg;

						//For N2
						dA = -dA;
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxI[3][k] = dsfg;
        	                        	
						//dGJI/dxi for N2
						dsfg = exp(-eta*ANN*ANN)*(-2*eta*ANN*dA*molec[mi].fcv[mj]);
						dgdxJ[3][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN + nsfg3point)	//point g3
				{
					//r =abs( molec[mi].neighdist[mj]);
					for (int k = 0; k < 3; k++)
					{
						dsfg = Rskappa*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj]*sin(Rskappa*molec[mi].neighdist[mj]) 
							+ cos(Rskappa*molec[mi].neighdist[mj])*molec[mi].dfcvdx[mj][k];
						dgdxI[0][k] = dsfg;
						dsfg =  dsfg;
						dgdxJ[0][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN + nsfg3point + nsfg3CO)	//CO g3
				{
						
					for (int k = 0; k < 3; k++)
					{
						CCO = Rskappa * rikjl/(magrik*magrjl);
						dC = Rskappa*(-nume1CO*new_molvec2CO[k] + nume2CO*new_molvec1CO[k])/denoCO;
			                        dsfg = -dC*sin(CCO)*molec[mi].fcv[mj] + cos(CCO)*molec[mi].dfcvdx[mj][k];
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
						dC = -dC;
						dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
						dgdxI[1][k] = dsfg;
						dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
						dgdxJ[1][k] = dsfg;
					}
				}
				else if (mk < nsfg2point + nsfg2CO + nsfg2NN + nsfg3point + nsfg3CO + nsfg3NN)	//NN g3
				{						
					for (int k = 0; k < 3; k++)
					{
						CNN = Rskappa * rmonp/(magrmo*magrnp);
						dC = Rskappa*(-nume1NN*new_molvec2NN[k] + nume2NN*new_molvec1NN[k])/denoNN;
	        	                        dsfg = cos(CNN)*molec[mi].dfcvdx[mj][k];
						dgdxI[0][k] = dsfg;
						dgdxJ[0][k] = dsfg;
	                	                dsfg = -dC*sin(CNN)*molec[mi].fcv[mj];
						dgdxI[2][k] = dsfg;
						dgdxJ[2][k] = dsfg;
				//For N2
						dC = -dC;
						dsfg = -dC*sin(CNN)*molec[mi].fcv[mj];
						dgdxI[3][k] = dsfg;
						dsfg = -dC*sin(CNN)*molec[mi].fcv[mj];
						dgdxJ[3][k] = dsfg;
					}
				}
				else
					cout << "ERROR, exceeded the limit" <<endl;
				for (int iex = 0; iex < parameter.nex;iex++)
				{
					for (int k = 0; k < 3; k++)
					{
						outi = outNN_index[iex]*parameter.nsfg + mk;
						qII = atoms[icenter].NNgrad[outi]/nmol;
						qJJ = atoms[jcenter].NNgrad[outi]/nmol;
						
						exvar[iex].dQ[atmi][k]+=qII*dgdxI[0][k]+qJJ*dgdxJ[0][k];
						exvar[iex].dQ[atmk][k]+=qII*dgdxI[1][k]+qJJ*dgdxJ[1][k];
						exvar[iex].dQ[atmm][k]+=qII*dgdxI[2][k]+qJJ*dgdxJ[2][k];
						exvar[iex].dQ[atmo][k]+=qII*dgdxI[3][k]+qJJ*dgdxJ[3][k];
					
					}//loop over xyz for calc dQ
				}//loop over extended variables
			}//loop over all symmetry functions
		}//loop over neighbor molecules
	}//loop over target molecules
}









//==================================================================================================================================================================
//stand alone only. here we calculate symmetry function in one function. 
void DAFED::symf_multi_vec_sa(Cmolpoint *molec, CParameter &parameter, int typeNum, int sfgNum){	//typeNum from symfType. sfgNum from i in nsfg loop
	int mi,mj, jmol;
        double Rskappa = parameter.RskappaLst[sfgNum];
	double eta = parameter.etaLst[sfgNum];
	double r;
	double dr;

	if (parameter.pointflag[sfgNum])	//point representation
	{
		if(parameter.g2g3flag[sfgNum])				//G2 symmetry function
		{
			double g2;
			for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
                		g2 = 0;
                		if (molec[mi].center==-1)exit(0);
                		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                		        jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index
	
	                        	dr = molec[mi].neighdist[mj]-Rskappa;
	                        	g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];
	                	}
	                	molec[mi].symf = g2;
	        	}
		}
		else					//G3 symmetry function
		{
			double g3;
			for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
				g3 = 0;
                		if (molec[mi].center==-1)exit(0);
                		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                        		jmol = molec[mi].neighbors[mj];

                       			r = abs(molec[mi].neighdist[mj]);
                        		g3 += cos(Rskappa*r)*molec[mi].fcv[mj];
                		}
                		molec[mi].symf = g3;
        		}
		}
	}
	else					//point-vector representation
	{
		double new_molvecTarget[3];		//cvecchunk1->pu_molvecall;
		double new_molvecNeighbour[3];
		double cosvec1;
		if(parameter.g2g3flag[sfgNum])					//G2 symmetry function
		{
			double g2;
        	        
			for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
				g2 = 0;                         //g2 start from 0 and adds up in the g2
                		if (molec[mi].center==-1)exit(1);	
				for (int k = 0; k < 3; k++)     //xyz component of the vector.
					new_molvecTarget[k] = molec[mi].vecpoint[typeNum][k];            //assigning the molvec of target here. typeNum should start from 0 and store xyz of vec in each		
				for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                	        	jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index
                        	
					for (int k = 0; k < 3; k++)             //xyz component of the vector for neighbor
                        		{
                        		        new_molvecNeighbour[k] = molec[jmol].vecpoint[typeNum][k];          //assigning the molvec
                        		}
                        		cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));          //cos(theta_ikjl)
                        		dr = cosvec1-Rskappa;                                                                     //dr is cosvec1-cosRs
                        		g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];                                                //sum of G2 for the target molecule
                		}
			
               		molec[mi].symf = g2;                                                                         //assigning sum of the g2 to each molecule 
			}
		}
		else						//G3 symmetry function
		{
			double g3;
			for(mi=0;mi<parameter.nmol;mi++){       //loop over all molecules
                		g3 = 0;                         //g2 start from 0 and adds up in the g2
                		if (molec[mi].center==-1)exit(1);	
				for (int k = 0; k < 3; k++)     //xyz component of the vector.
					new_molvecTarget[k] = molec[mi].vecpoint[typeNum][k];            //assigning the molvec of target here. typeNum should start from 0 and store xyz of vec in each		
                		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                        		jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index

                        		for (int k = 0; k < 3; k++)             //xyz component of the vector for neighbor
                        		{
                        		        new_molvecNeighbour[k] = molec[jmol].vecpoint[typeNum][k];          //assigning the molvec
                        		}
					cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));
					g3 += cos(Rskappa*cosvec1)*molec[mi].fcv[mj];
                		}
                	molec[mi].symf = g3;                                                           //assigning sum of the g2 to each molecule 
			}
		}
	}
}
//stand alone only. short version of one above. Working.  
void DAFED::symf_multi_vec_sa_short_test(CAtom *atoms, Cmolpoint *molec, CParameter &parameter){	//typeNum from symfType. sfgNum from i in nsfg loop
	unsigned short int mi,mj, jmol,i;
        double Rskappa, eta;
	double dr;
	double G;
	double new_molvecTarget[3] ={0,0,0};		//cvecchunk1->pu_molvecall;
	double new_molvecNeighbour[3]={0,0,0};
	double cosvec1;
	double neidist;
	double fcv;
	short typeNum;
	short old_typeNum = -1;
	short nsfg=parameter.nsfg;
	unsigned short int nmol = parameter.nmol;
	int n_nei = -1;
	for(mi=0;mi<nmol;++mi)       //loop over all molecules
	{
		n_nei = molec[mi].n_neighbors;	
		for (mj = 0; mj < n_nei;++mj)	//loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
		{
			jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index
			neidist = molec[mi].neighdist[mj];
			fcv = molec[mi].fcv[mj];

			for (i = 0; i < nsfg; ++i)
			{
        			Rskappa = parameter.RskappaLst[i];
				eta = parameter.etaLst[i];
				G = 0;
				typeNum = parameter.symftypeLst[i];
				
				if (parameter.pointflag[i])	//point representation
				{
					if(parameter.g2g3flag[i])				//G2 symmetry function
					{
						dr = neidist-Rskappa;
						G = exp(-eta*dr*dr)*fcv;
					}
					else							//G3 symf
					{
                		       		dr = abs(neidist);
                        			G = cos(Rskappa*dr)*fcv;
                			}
				}
				else					//point vector
				{
					if(old_typeNum != typeNum)
					{
					for (short k = 0; k < 3; ++k)     //xyz component of the vector.
					{	
						new_molvecTarget[k] = molec[mi].vecpoint[typeNum][k];	//assigning the molvec of target here. typeNum should start from 0 and store xyz of vec in each
                        			new_molvecNeighbour[k] = molec[jmol].vecpoint[typeNum][k];	//assigning the molvec
					}
					old_typeNum = typeNum;
					}
					
					if(parameter.g2g3flag[i])					//G2 symmetry function
					{				
                        			cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));          //cos(theta_ikjl)
                        			dr = cosvec1-Rskappa;                                                                     //dr is cosvec1-cosRs
                        			G = exp(-eta*dr*dr)*fcv;                                                //sum of G2 for the target molecule
					}
					else						//G3 symmetry function
					{
					
						cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));
						G = cos(Rskappa*cosvec1)*fcv;
					}
	                	}
				atoms[molec[mi].center].sfg[i] += G;
	        	}
		}
	}
}


//stand alone only. calculate derivative of symf.
void DAFED::get_total_derivatives_molvec_vec_sa(CAtom *atoms, Cmolpoint *molec, int &natoms,CParameter &parameter, EXvariables *exvar){
	int mi,mj,mk;
	int jmol;
        double Rskappa;
        double eta;
	double dsfg;
	double dr;
	double new_molvecTarget[parameter.symftype][3];
	double new_molvecNeighbour[parameter.symftype][3];
	double CCO,dC, ACO, dA;
        double rikjl, magrik, magrjl, nume1,nume2, deno;
	int nsfg = parameter.nsfg;
	double qII,qJJ;
	int outi;
	int outNN_index[2];
	double dgdxICenter[3],dgdxI_1[3],dgdxI_2[3],dgdxJCenter[3],dgdxJ_1[3],dgdxJ_2[3];
	double nmol = parameter.nmol;
	int typeNum = -1;
	outNN_index[0] =0;		//output of urea I in the NN
	//outNN_index[1] =2;		//output of urealiq in the NN
	outNN_index[1] =1;		//output of urea IV in the NN
	for(int iex=0;iex<parameter.nex;iex++){
	        for(int ti=0;ti<natoms;ti++){	//	******This code does not work on parallel because of this line. If you want to run in parallel, change this to local + ghost ***********
			for(int i=0;i<3;i++){
				exvar[iex].dQ[ti][i] = 0.0;
			}
		}
	}
	for(mi=0;mi<nmol;mi++)		//loop over all molecules original was parameter.nmol, but changed to 10 for testing.
	{
		if (molec[mi].center ==-1)break;
		
		int icenter = molec[mi].center;
		
		for (int sfg = 2; sfg<parameter.symftype;sfg++)
		{
			typeNum = parameter.symftypeLst[sfg];
			for (int k = 0; k < 3; k++)		//get molvector of target molecule
			{             
	                	new_molvecTarget[sfg][k] = atoms[icenter].vecpoint[typeNum][k];
			}
		}
		for (mj = 0; mj < molec[mi].n_neighbors;mj++)		//loop over all the neighbor molecules.
		{ 
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			int jcenter = molec[jmol].center;                 //center of the molec jmol id. 
			if (icenter==jcenter)                     //Skip if target and neighbor is the same molecule
				continue;	
			for (int sfg = 2; sfg<parameter.symftype;sfg++)
                	{
				typeNum = parameter.symftypeLst[sfg];
				for (int k = 0; k < 3; k++)
				{
					new_molvecNeighbour[sfg][k] = atoms[jcenter].vecpoint[typeNum][k];
				}
			}
			for (mk = 0; mk < nsfg; mk++)
			{
				
				

				for(int k = 0;k<3;k++)  //for x,y,z
				{
				dgdxICenter[k] = 0;
				dgdxI_1[k] = 0;
				dgdxI_2[k] = 0;
				dgdxJCenter[k] = 0;
				dgdxJ_1[k] = 0;
				dgdxJ_2[k] = 0;
				}
				Rskappa = parameter.RskappaLst[mk];
	                        eta = parameter.etaLst[mk];
			
			
				if (parameter.pointflag[mk])	//point representation
				{
					if(parameter.g2g3flag[mk])				//G2 symmetry function
					{
						dr = molec[mi].neighdist[mj]-Rskappa;
	                                        for (int k = 0; k < 3; k++)
		                                {
							dsfg = exp(-eta*dr*dr)*(2.0*eta*dr*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);
                         	                	dgdxICenter[k] = dsfg;//0 because this is for the point
                                	        	dgdxJCenter[k] = dsfg;//0 because this is for the point
                                        	}
					}
					else
					{
						for (int k = 0; k < 3; k++)
						{
							dsfg = Rskappa*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj]*sin(Rskappa*molec[mi].neighdist[mj]) 
								+ cos(Rskappa*molec[mi].neighdist[mj])*molec[mi].dfcvdx[mj][k];
							dgdxICenter[k] = dsfg;//0 because this is for the point
							dgdxJCenter[k] = dsfg;//0 because this is for the point


						}
					}
				}
				else
				{
					rikjl = dot3(new_molvecTarget[mk],new_molvecNeighbour[mk]);
					magrik = (len3(new_molvecTarget[mk]));
					magrjl = (len3(new_molvecNeighbour[mk]));
					nume1 = magrik*magrjl;
					nume2 = rikjl*magrjl/magrik;
					deno = magrik*magrjl*magrik*magrjl;

					if(parameter.g2g3flag[mk])				//G2 symmetry function
					{
						if (parameter.vectype[parameter.symftypeLst[mk]*2] == parameter.vectype[0])	//if element used for vec is same as center...
							for (int k = 0; k < 3; k++)
							{
								ACO = rikjl/(magrik*magrjl) - Rskappa;

								dA = (-nume1*new_molvecNeighbour[mk][k] + nume2*new_molvecTarget[mk][k])/deno;				//derivative of A
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);           //derivative of g by xi,yi,zi
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
	
								dA = -dA;					//When we have a derivative of G2I or G2J with atmk, dA will be just flipped.
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}

						else if (parameter.vectype[parameter.symftypeLst[mk]*2 + 1] == parameter.vectype[0])//if another element used for vec is same as center...flipped dsfg
						{
							for (int k = 0; k < 3; k++)
							{
								ACO = rikjl/(magrik*magrjl) - Rskappa;
	
								dA = (-nume1*new_molvecNeighbour[mk][k] + nume2*new_molvecTarget[mk][k])/deno;				//derivative of A
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
	
								dA = -dA;					//When we have a derivative of G2I or G2J with atmk, dA will be just flipped.
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);           //derivative of g by xi,yi,zi
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}	
						
						}

						else	//if vector elements are not used as center. 
						{
							for (int k = 0; k < 3; k++)
							{
								ACO = rikjl/(magrik*magrjl) - Rskappa;

								dA = (-nume1*new_molvecNeighbour[mk][k] + nume2*new_molvecTarget[mk][k])/deno;				//derivative of A
								dsfg = exp(-eta*ACO*ACO)*(molec[mi].dfcvdx[mj][k]);
								dgdxICenter[k] = dsfg;		//this is for the center
								dgdxJCenter[k] = dsfg;
		
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);
								dgdxI_1[k] = dsfg;		//this is for the vec ele 1
								dgdxJ_1[k] = dsfg;
		
								dA = -dA;					//When we have a derivative of G2I or G2J with atmk, dA will be just flipped.
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);           //derivative of g by xi,yi,zi
								dgdxI_2[k] = dsfg;		//this is for the vec ele 2
								dgdxJ_2[k] = dsfg;
							}							
						
						}
					}
					else
					{
						if (parameter.vectype[parameter.symftypeLst[mk]*2] == parameter.vectype[0])	//if element used for vec is same as center...
						{
							for (int k = 0; k < 3; k++)
							{
								CCO = Rskappa * rikjl/(magrik*magrjl);
								dC = Rskappa*(-nume1*new_molvecNeighbour[mk][k] + nume2*new_molvecTarget[mk][k])/deno;
			                        		dsfg = -dC*sin(CCO)*molec[mi].fcv[mj] + cos(CCO)*molec[mi].dfcvdx[mj][k];
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
								dC = -dC;
								dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}

						}
						else if (parameter.vectype[parameter.symftypeLst[mk]*2 + 1] == parameter.vectype[0])//if another element used for vec is same as center...flipped dsfg
						{
							for (int k = 0; k < 3; k++)
							{
								CCO = Rskappa * rikjl/(magrik*magrjl);
								dC = Rskappa*(-nume1*new_molvecNeighbour[mk][k] + nume2*new_molvecTarget[mk][k])/deno;
								dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
								dC = -dC;
			                	        	dsfg = -dC*sin(CCO)*molec[mi].fcv[mj] + cos(CCO)*molec[mi].dfcvdx[mj][k];
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}
						}
						else
						{
							for (int k = 0; k < 3; k++)
							{
								CCO = Rskappa * rikjl/(magrik*magrjl);
								dC = Rskappa*(-nume1*new_molvecNeighbour[mk][k] + nume2*new_molvecTarget[mk][k])/deno;
		        	                        	dsfg = cos(CCO)*molec[mi].dfcvdx[mj][k];
								dgdxICenter[k] = dsfg;
								dgdxJCenter[k] = dsfg;
	                		                	dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
								//For N2
								dC = -dC;
								dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}
						}
					}	
				}
				
				for (int iex = 0; iex < parameter.nex;iex++)
				{
					for (int k = 0; k < 3; k++)
					{
						outi = outNN_index[iex]*parameter.nsfg + mk;
						qII = atoms[icenter].NNgrad[outi]/nmol;
						qJJ = atoms[jcenter].NNgrad[outi]/nmol;
						
						exvar[iex].dQ[0][k]+=qII*dgdxICenter[k]+qJJ*dgdxJCenter[k];//center
						exvar[iex].dQ[parameter.vectype[mk]*2][k]+=qII*dgdxI_1[k]+qJJ*dgdxJ_1[k];//center
						exvar[iex].dQ[parameter.vectype[mk]*2+1][k]+=qII*dgdxI_2[k]+qJJ*dgdxJ_2[k];//center

					}//loop over xyz for calc dQ
				}//loop over extended variables
			}//loop over all symmetry functions
		}//loop over neighbor molecules
	}//loop over target molecules
}







//from here lammps only 
//Need to make symf for lmp in here
void DAFED::symf_multi_vec_lmp(CAtom *molecules, Cmolpoint *molec, CParameter &parameter, CNeighvariables &lmpneigh, int typeNum, int sfgNum){	//typeNum from symfType. sfgNum from i in nsfg loop
	int mi,mj, jmol;
        double Rskappa = parameter.RskappaLst[sfgNum];
	double eta = parameter.etaLst[sfgNum];
	double r;
	double dr;

	if (parameter.pointflag[sfgNum])	//point representation
	{
		if(parameter.g2g3flag[sfgNum])				//G2 symmetry function
		{
			double g2;
			for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
                		g2 = 0;
                		if (molec[mi].center==-1)exit(0);
                		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                		        jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index
	
	                        	dr = molec[mi].neighdist[mj]-Rskappa;
	                        	g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];
	                	}
	                	molec[mi].symf = g2;
	        	}
		}
		else					//G3 symmetry function
		{
			double g3;
			for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
				g3 = 0;
                		if (molec[mi].center==-1)exit(0);
                		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                        		jmol = molec[mi].neighbors[mj];

                       			r = abs(molec[mi].neighdist[mj]);
                        		g3 += cos(Rskappa*r)*molec[mi].fcv[mj];
                		}
                		molec[mi].symf = g3;
        		}
		
		}
		
	}
	else					//point-vector representation
	{
		double new_molvecTarget[3];		//cvecchunk1->pu_molvecall;
		double new_molvecNeighbour[3];
		double cosvec1;
		if(parameter.g2g3flag[sfgNum])					//G2 symmetry function
		{
			double g2;
			for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
				g2 = 0;                         //g2 start from 0 and adds up in the g2
                		if (molec[mi].center==-1)exit(1);	
				for (int k = 0; k < 3; k++)     //xyz component of the vector.
					new_molvecTarget[k] = molecules[molec[mi].center].vecpoint[typeNum][k];		//assigning the molvec of target here. typeNum should start from 0 and store xyz of vec in each		
				for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                	        	jmol = molec[mi].neighbors[mj];         //assigning the neighbor atom's id
                        	
					for (int k = 0; k < 3; k++)             //xyz component of the vector for neighbor
                        		{
                        		        new_molvecNeighbour[k] = molecules[jmol].vecpoint[typeNum][k];          //assigning the molvec
                        		}
                        		cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));          //cos(theta_ikjl)
                        		dr = cosvec1-Rskappa;                                                                     //dr is cosvec1-cosRs
                        		g2 += exp(-eta*dr*dr)*molec[mi].fcv[mj];                                                //sum of G2 for the target molecule

                		}
			
               		molec[mi].symf = g2;                                                                         //assigning sum of the g2 to each molecule 
			}
		}
		else						//G3 symmetry function
		{
			double g3;
			for(mi=0;mi<lmpneigh.lmpmol;mi++){       //loop over all molecules
                		g3 = 0;                         //g2 start from 0 and adds up in the g2
                		if (molec[mi].center==-1)exit(1);	
				for (int k = 0; k < 3; k++)     //xyz component of the vector.
					new_molvecTarget[k] = molecules[molec[mi].center].vecpoint[typeNum][k];            //assigning the molvec of target here. typeNum should start from 0 and store xyz of vec in each		
                		for (mj = 0; mj < molec[mi].n_neighbors;mj++){  //loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
                        		jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index

                        		for (int k = 0; k < 3; k++)             //xyz component of the vector for neighbor
                        		{
                        		        new_molvecNeighbour[k] = molecules[jmol].vecpoint[typeNum][k];          //assigning the molvec
                        		}
					cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));
					g3 += cos(Rskappa*cosvec1)*molec[mi].fcv[mj];
                		}
                	molec[mi].symf = g3;                                                           //assigning sum of the g2 to each molecule 
			}
		}
	}
}

//this is the latest version of lammps symf calculation
void DAFED::symf_multi_vec_lmp_short_test(CAtom *atoms, Cmolpoint *molec, CParameter &parameter, CNeighvariables &lmpneigh){
	unsigned short int mi,mj, jmol,i;
        double Rskappa, eta;
	double dr;
	double G;
	double new_molvecTarget[3]={0,0,0};		//cvecchunk1->pu_molvecall;
	double new_molvecNeighbour[3]={0,0,0};
	double cosvec1;
	double neidist;
	double fcv;
	short typeNum;
	short old_typeNum = -1;
	short nsfg=parameter.nsfg;
	int n_nei = -1;
	for(mi=0;mi<lmpneigh.lmpmol;++mi)       //loop over all molecules
	{
		n_nei = molec[mi].n_neighbors;	
		for (mj = 0; mj < n_nei;++mj)	//loop over all the neighbor molecules.      This can be improved if I use the molec.neighbors
		{
			jmol = molec[mi].neighbors[mj];         //assigning the neighbor molecule index
			neidist = molec[mi].neighdist[mj];
			fcv = molec[mi].fcv[mj];
                        
			for (i = 0; i < nsfg; ++i)
			{
        			Rskappa = parameter.RskappaLst[i];
				eta = parameter.etaLst[i];
				G = 0;
				typeNum = parameter.symftypeLst[i];
				
				if (parameter.pointflag[i])	//point representation
				{
					if(parameter.g2g3flag[i])				//G2 symmetry function
					{
						dr = neidist-Rskappa;
						G = exp(-eta*dr*dr)*fcv;
					}
					else							//G3 symf
					{
                		       		dr = abs(neidist);
                        			G = cos(Rskappa*dr)*fcv;
                			}
				}
				else					//point vector
				{
					if(old_typeNum != typeNum)
					{
					for (short k = 0; k < 3; ++k)     //xyz component of the vector.
					{	
						new_molvecTarget[k] = atoms[molec[mi].center].vecpoint[typeNum][k];	//assigning the molvec of target here. typeNum should start from 0 and store xyz of vec in each
                        			new_molvecNeighbour[k] = atoms[jmol].vecpoint[typeNum][k];;	//assigning the molvec
					}
					old_typeNum = typeNum;
					}
					if(parameter.g2g3flag[i])					//G2 symmetry function
					{				
                        			cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));          //cos(theta_ikjl)
                        			dr = cosvec1-Rskappa;                                                                     //dr is cosvec1-cosRs
                        			G = exp(-eta*dr*dr)*fcv;                                                //sum of G2 for the target molecule
					}
					else						//G3 symmetry function
					{
						cosvec1 = dot3(new_molvecTarget,new_molvecNeighbour)/(len3(new_molvecTarget)*len3(new_molvecNeighbour));
						G = cos(Rskappa*cosvec1)*fcv;	
					}

	                	}
				atoms[molec[mi].center].sfg[i] += G;
	        	}
		}
	}
}



//This is the latest version of how to calculate derivative of each symmetry function.
void DAFED::get_total_derivatives_molvec_vec_lmp(CAtom *atoms, Cmolpoint *molec, int &natoms,CParameter &parameter, EXvariables *exvar,CNeighvariables &lmpneigh){
	int mi,mj,mk;
	int jmol;
        double Rskappa;
        double eta;
	double dsfg;
	double dr;
	double new_molvecTarget[parameter.symftype][3];
	double new_molvecNeighbour[parameter.symftype][3];
	double CCO,dC, ACO, dA;
        double rikjl, magrik, magrjl, nume1,nume2, deno;
	int nsfg = parameter.nsfg;
	double qII,qJJ;
	int outi;
	int outNN_index[2];
	double dgdxICenter[3],dgdxI_1[3],dgdxI_2[3],dgdxJCenter[3],dgdxJ_1[3],dgdxJ_2[3];
	int typeNum = -1;
	outNN_index[0] =0;		//output of urea I in the NN
	//outNN_index[1] =2;		//output of urealiq in the NN
	outNN_index[1] =1;		//output of urea IV in the NN
	double nmol = parameter.nmol;
	double lmpmol = lmpneigh.lmpmol;
	for(int iex=0;iex<parameter.nex;iex++){
	        for(int ti=0;ti<lmpneigh.nall;ti++){
			for(int i=0;i<3;i++){
				exvar[iex].dQ[ti][i] = 0.0;
			}
		}
	}
	for(mi=0;mi<lmpmol;mi++)		//loop over all molecules original was parameter.nmol, but changed to 10 for testing.
	{
		if (molec[mi].center ==-1)break;
		
		int icenter = molec[mi].center;
		
		for (int sfg = 1; sfg<parameter.symftype;sfg++)
		{
			//typeNum = parameter.symftypeLst[sfg];
			for (int k = 0; k < 3; k++)		//get molvector of target molecule
			{             
	                	new_molvecTarget[sfg][k] = atoms[icenter].vecpoint[sfg][k];
			}
		}
		for (mj = 0; mj < molec[mi].n_neighbors;mj++)		//loop over all the neighbor molecules.
		{ 
			jmol = molec[mi].neighbors[mj];		//assigning the neighbor molecule index
			int jcenter = jmol;                 //center of the molec jmol id. 
			if (icenter==jcenter)                     //Skip if target and neighbor is the same molecule
				continue;	
			for (int sfg = 1; sfg<parameter.symftype;sfg++)
                	{
				for (int k = 0; k < 3; k++)
				{
					new_molvecNeighbour[sfg][k] = atoms[jcenter].vecpoint[sfg][k];
				}
			}

			for (mk = 0; mk < nsfg; mk++)
			{
				for(int k = 0;k<3;k++)  //for x,y,z
				{
				dgdxICenter[k] = 0;	//center atom
				dgdxI_1[k] = 0;		//first atom of the vector
				dgdxI_2[k] = 0;		//second atom of the vector
				dgdxJCenter[k] = 0;	//same but for neighbor
				dgdxJ_1[k] = 0;
				dgdxJ_2[k] = 0;
				}
				Rskappa = parameter.RskappaLst[mk];
	                        eta = parameter.etaLst[mk];
			
			
				if (parameter.pointflag[mk])	//point representation
				{
					typeNum = -1;	//used in dQ calc. if it is point, it should not consider neighbors
					if(parameter.g2g3flag[mk])				//G2 symmetry function
					{
						dr = molec[mi].neighdist[mj]-Rskappa;
	                                        for (int k = 0; k < 3; k++)
		                                {
							dsfg = exp(-eta*dr*dr)*(2.0*eta*dr*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);
                         	                	dgdxICenter[k] = dsfg;//0 because this is for the point
                                	        	dgdxJCenter[k] = dsfg;//0 because this is for the point
                                        	}
					}
					else
					{
						for (int k = 0; k < 3; k++)
						{
							dsfg = Rskappa*molec[mi].diff[mj][k]/molec[mi].neighdist[mj]*molec[mi].fcv[mj]*sin(Rskappa*molec[mi].neighdist[mj]) 
								+ cos(Rskappa*molec[mi].neighdist[mj])*molec[mi].dfcvdx[mj][k];
							dgdxICenter[k] = dsfg;//0 because this is for the point
							dgdxJCenter[k] = dsfg;//0 because this is for the point


						}
					}
				}
				else
				{
					typeNum = parameter.symftypeLst[mk];
					rikjl = dot3(new_molvecTarget[typeNum],new_molvecNeighbour[typeNum]);
					magrik = (len3(new_molvecTarget[typeNum]));
					magrjl = (len3(new_molvecNeighbour[typeNum]));
					nume1 = magrik*magrjl;
					nume2 = rikjl*magrjl/magrik;
					deno = magrik*magrjl*magrik*magrjl;

			/////Maybe I can have a variable call vectypeEle1 = parameter.vectype[mk]*2;
			//				      vectypeEle2 = parameter.vectype[mk]*2+1
					if(parameter.g2g3flag[mk])				//G2 symmetry function
					{
						if (parameter.vectype[parameter.symftypeLst[mk]*2] == parameter.vectype[0])	//if element used for vec is same as center...
							for (int k = 0; k < 3; k++)
							{
								ACO = rikjl/(magrik*magrjl) - Rskappa;

								dA = (-nume1*new_molvecNeighbour[typeNum][k] + nume2*new_molvecTarget[typeNum][k])/deno;				//derivative of A
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);           //derivative of g by xi,yi,zi
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
	
								dA = -dA;					//When we have a derivative of G2I or G2J with atmk, dA will be just flipped.
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}

						else if (parameter.vectype[parameter.symftypeLst[mk]*2 + 1] == parameter.vectype[0])//if another element used for vec is same as center...flipped dsfg
						{
							for (int k = 0; k < 3; k++)
							{
								ACO = rikjl/(magrik*magrjl) - Rskappa;
	
								dA = (-nume1*new_molvecNeighbour[typeNum][k] + nume2*new_molvecTarget[typeNum][k])/deno;				//derivative of A
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
	
								dA = -dA;					//When we have a derivative of G2I or G2J with atmk, dA will be just flipped.
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj] + molec[mi].dfcvdx[mj][k]);           //derivative of g by xi,yi,zi
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}	
						
						}

						else	//if vector elements are not used as center. 
						{
							for (int k = 0; k < 3; k++)
							{
								ACO = rikjl/(magrik*magrjl) - Rskappa;

								dA = (-nume1*new_molvecNeighbour[typeNum][k] + nume2*new_molvecTarget[typeNum][k])/deno;				//derivative of A
								dsfg = exp(-eta*ACO*ACO)*(molec[mi].dfcvdx[mj][k]);
								dgdxICenter[k] = dsfg;		//this is for the center
								dgdxJCenter[k] = dsfg;
		
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);
								dgdxI_1[k] = dsfg;		//this is for the vec ele 1
								dgdxJ_1[k] = dsfg;
		
								dA = -dA;					//When we have a derivative of G2I or G2J with atmk, dA will be just flipped.
								dsfg = exp(-eta*ACO*ACO)*(-2*eta*ACO*dA*molec[mi].fcv[mj]);           //derivative of g by xi,yi,zi
								dgdxI_2[k] = dsfg;		//this is for the vec ele 2
								dgdxJ_2[k] = dsfg;
							}							
						
						}
					}
					else	//vector representation
					{
						if (parameter.vectype[parameter.symftypeLst[mk]*2] == parameter.vectype[0])	//if element used for vec is same as center...
						{
							for (int k = 0; k < 3; k++)
							{
								CCO = Rskappa * rikjl/(magrik*magrjl);
								dC = Rskappa*(-nume1*new_molvecNeighbour[typeNum][k] + nume2*new_molvecTarget[typeNum][k])/deno;
			                        		dsfg = -dC*sin(CCO)*molec[mi].fcv[mj] + cos(CCO)*molec[mi].dfcvdx[mj][k];
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
								dC = -dC;
								dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}

						}
						else if (parameter.vectype[parameter.symftypeLst[mk]*2 + 1] == parameter.vectype[0])//if another element used for vec is same as center...flipped dsfg
						{
							for (int k = 0; k < 3; k++)
							{
								CCO = Rskappa * rikjl/(magrik*magrjl);
								dC = Rskappa*(-nume1*new_molvecNeighbour[typeNum][k] + nume2*new_molvecTarget[typeNum][k])/deno;
								dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
								dC = -dC;
			                	        	dsfg = -dC*sin(CCO)*molec[mi].fcv[mj] + cos(CCO)*molec[mi].dfcvdx[mj][k];
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}
						}
						else	//if neither of vec elements are center
						{
							for (int k = 0; k < 3; k++)
							{
								CCO = Rskappa * rikjl/(magrik*magrjl);
								dC = Rskappa*(-nume1*new_molvecNeighbour[typeNum][k] + nume2*new_molvecTarget[typeNum][k])/deno;
		        	                        	dsfg = cos(CCO)*molec[mi].dfcvdx[mj][k];
								dgdxICenter[k] = dsfg;
								dgdxJCenter[k] = dsfg;
	                		                	dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_1[k] = dsfg;
								dgdxJ_1[k] = dsfg;
								dC = -dC;
								dsfg = -dC*sin(CCO)*molec[mi].fcv[mj];
								dgdxI_2[k] = dsfg;
								dgdxJ_2[k] = dsfg;
							}//end of xyz
						}//end of center and 2 vectors are different case
					}	//end of G3 point vector ysmmetry function
				}		//end of 1 symmetry function or I can say end of point vector symmetry function
				for (int iex = 0; iex < parameter.nex;iex++)
				{
					for (int k = 0; k < 3; k++)
					{
						outi = outNN_index[iex]*parameter.nsfg + mk;
						qII = atoms[icenter].NNgrad[outi]/nmol;
						qJJ = atoms[jcenter].NNgrad[outi]/nmol;
						
						exvar[iex].dQ[icenter][k]+=qII*dgdxICenter[k]+qJJ*dgdxJCenter[k];//center
						
						
						if(typeNum!=-1)
						{
						exvar[iex].dQ[molec[mi].vector_ID[typeNum*2]][k]+=qII*dgdxI_1[k]+qJJ*dgdxJ_1[k];//first atom of the vec
						exvar[iex].dQ[molec[mi].vector_ID[typeNum*2+1]][k]+=qII*dgdxI_2[k]+qJJ*dgdxJ_2[k];//second atom of the vec 
						}
						
						for(int ir=0;ir<3;ir++)
						{        //loop over x,y,z for position, add up virial
							exvar[iex].dvirial[ir][k] += (qII*(dgdxICenter[k]+dgdxI_1[k]+dgdxI_2[k])
								+qJJ*(dgdxJCenter[k]+dgdxJ_1[k]+dgdxJ_2[k]))*molec[mi].diff[mj][ir]*0.5;
                                        	}//loop over xyz for virial
					}//loop over xyz for calc dQ
				}//loop over extended variables
			}//loop over all symmetry functions
		//cout <<"End of each neigh loop " << endl;
		}//loop over neighbor molecules
		//cout <<"End of each target loop " << endl;
	}//loop over target molecules
}
