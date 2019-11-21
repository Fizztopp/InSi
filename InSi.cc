/**
 *	TIGHT-BINDING MODEL FOR ONE_DIMENSIONAL INDIUM WIRES (InSi)
 *  Copyright (C) 2019, Gabriel E. Topp
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 *  02111-1307, USA.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <assert.h>
#include <iterator>
#include <sstream>
#include <string>
#include <algorithm>


// PARAMETERS ##########################################################

#define PI        3.14159265359

// Intrinsic parameters
// Electronic
#define NATOM     8                           						    // Number of orbitals (atoms) in unit cell
#define NN        1024                                                  // # of k-points
#define eO        0.056                                                 // Onsite energy 1
#define eI        0.089                                                 // Onsite energy 2
#define tO        -0.290                                                // Hopping 1            above Tc: -0.419  below Tc: -0.350 (ROTARY)
#define tOd       -0.545                                                // Hopping 1 dashed     above Tc: -0.419  below Tc: -0.501 (ROTARY)
#define tI1       0.100                    					            // Hopping 2            above Tc: 0.284   below Tc: 0.145  (SHEAR)
#define tI1d      0.550                                                 // Hopping 2 dashed     above Tc: 0.284   below Tc: 0.481  (SHEAR)
#define tI2       -0.104												// hopping 3
#define tIO       0.147													// hopping 4
#define BETA      250.                       					     	// inverse temperature  300 K -> 40, 40 K --> 250

// Numerical paramters
#define mu_init   -0.1											     	// initially guessed chemical potential
#define dev       1e-8                    					        	// exit deviation for while loop in groundstate() (m)
#define DELTA     1e-5												    // correction prefactor for chemical potential in groundstate(

// Peierls driving
#define w_peierls      0.188                                            // Frequency of Applied Field (in eV)
#define Ax_peierls     0.1                                              // Amplitude of Applied Field in x-direction
#define Ay_peierls     0.0                                              // Amplitude of Applied Field in y-direction
#define SIGMA          274.0                                            // 
#define DELAY          1520.0                                         	// Mean value of Gauss

// Propagation parameters
#define starttime 0.0 
#define endtime   3040.0
#define TIMESTEPS 1e4                                                   

#ifndef NO_MPI
    #include <mpi.h>
#endif

using namespace std;

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);

// DEFINITION OF FUNCTIONS #############################################

//LAPACK (Fortran 90) functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//routine to find eigensystem of matrix
extern "C" {
/** 
 *  Computes the eigenvalues and, optionally, the eigenvectors for a Hermitian matrices H
 */
    void zheev_(char* jobz, char* uplo, int* N, cdouble* H, int* LDA, double* W, cdouble* work, int* lwork, double* rwork, int *info);
}
//'N','V':  Compute eigenvalues only, and eigenvectors
char    jobz = 'V';       
//'U','L':  Upper, Lower triangle of H is stored 
char    uplo = 'U';  
// The order of the matrix H.  NATOM >= 0
int     matsize = NATOM;    
// The leading dimension of the array H.  lda >= max(1, NATOM)
int     lda = NATOM;             
// The length of the array work.  lwork  >= max(1,2* NATOM-1)
int     lwork = 2*NATOM-1;    
// dimension (max(1, 3* NATOM-2))
double  rwork[3*NATOM-2];  
// dimension (MAX(1,LWORK))
cdouble work[2*NATOM-1];  
// Info
int	    info;


void diagonalize(cvec &Hk, dvec &evals)
{
/**
 *  Diagonalization of matrix Hk. Stores eigenvalues in real vector evals and eigenvectors in complex vector Hk
 *  -Hk: Complex vector[NATOM x NATOM] to store Hamiltonian --> transformation matrices
 * 	-evals: Real vector[NATOM] to store eigenvalues
 */
    zheev_(&jobz, &uplo, &matsize, &Hk[0], &lda, &evals[0], &work[0], &lwork, &rwork[0], &info);
	assert(!info);
}


//INLINE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline int fq(int i, int j, int N)
/**
 *  MAT[i,j] = Vec[fq(i,j,N)] with row index i and column index j
 */
{
    return i*N+j;
}


inline double delta(int a, int b)
/**
 *  Delta function
 */
{
	if (a==b)
		return 1.;
	else
		return 0.;
}


inline double Ax_t(double time)
{
/**
 *	Peierls field for electrons in x-direction:
 *  -time: Real time coordinate
 */
    return Ax_peierls*sin(w_peierls*time)*exp(-0.5*pow((time-DELAY)/SIGMA,2.));
}


inline double Ay_t(double time)
{
/**
 *	Peierls field for electrons in y-direction:
 *  -time: Real time coordinate
 */
    return Ay_peierls*sin(w_peierls*time)*exp(-0.5*pow((time-DELAY)/SIGMA,2.));
}


inline double fermi(double energy, double mu)
{
/**
 *	Fermi distribution:
 *	-energy: Energy eigenvalue
 *	-mu: Chemical potential
 */
    return 1./(exp((energy-mu)*BETA) + 1.);
}


inline double gauss(double time, double delay, double sigma)
/**
 *	Normalized Gauss distribution
 *	-time: time coordinate
 *	-delay: mean expectation value
 *	-sigma: standard deviation 
 **/
{
	return 1./(sigma*sqrt(2.*PI))*exp(-0.5*pow((time-delay)/sigma,2.));
}


// VOID FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template <class Vec>
void times(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product of quadratic matrices: $C = A \cdot B$
 */
{
    int dim = sqrt(A.size());
	Vec TEMP(dim*dim);
    // Transposition gives speed up due to avoided line break
	for(int i=0; i<dim; i++) {
	    for(int j=0; j<dim; j++) {
		    TEMP[fq(j,i,dim)] = B[fq(i,j,dim)];
		   }
    }
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
				C[fq(i,j,dim)] += A[fq(i,k,dim)]*TEMP[fq(j,k,dim)]; 
			}
		}
	}	
}


template <class Vec>
void times_dn(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product with Hermitian conjugation of first factor: $C = A^\dagger \cdot B$
 */
{
	int dim = sqrt(A.size());
	Vec TEMP1(dim*dim);
	Vec TEMP2(dim*dim);
	// Transposition gives speed up due to avoided line break
	for(int i=0; i<dim; i++) {
		for(int j=0; j<dim; j++) {
			TEMP1[fq(j,i,dim)] = A[fq(i,j,dim)];
			TEMP2[fq(j,i,dim)] = B[fq(i,j,dim)];
		}
	}		
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
				C[fq(i,j,dim)] += conj(TEMP1[fq(i,k,dim)])*TEMP2[fq(j,k,dim)];
			}
		}
	}		
}


template <class Vec>
void times_nd(Vec &A, Vec &B, Vec &C)
/**
 *	Matrix product with Hermitian conjugation of second factor: $C = A \cdot B^\dagger$
 */
{
	int dim = sqrt(A.size());	
	for(int i=0; i<dim; ++i)
	{
		for(int j=0; j<dim; ++j)
		{
			C[fq(i,j,dim)] = 0.;
			for(int k=0; k<dim; ++k)
			{
					C[fq(i,j,dim)] += A[fq(i,k,dim)]*conj(B[fq(j,k,dim)]);
			}
		}
	}	
}


void set_Hk(double k, cvec &Hk, double time)
/**
 *	Set time-dependent Hamiltonian matrix with Peierls field
 *  -k: Real vector of the reciprocal space
 *  -Hk: Complex vector[NATOM*NATOM] to store Hamiltonian
 *  -time: time variable
 */
{
	for(int i=0; i<8; ++i)
            for(int j=0; j<8; ++j)
                Hk[fq(i,j,8)] = 0.;
    
    // Diagonal elements    
	Hk[fq(0,0,8)] = eO;	
	Hk[fq(1,1,8)] = eI; 
	Hk[fq(2,2,8)] = eI;
	Hk[fq(3,3,8)] = eO;
	Hk[fq(4,4,8)] = eO;
	Hk[fq(5,5,8)] = eI;
	Hk[fq(6,6,8)] = eI;
	Hk[fq(7,7,8)] = eO;	 
	
	// Hopping elements
	// Lower triagonal
	Hk[fq(7,0,8)] += -tOd*exp(-II*Ax_t(time))*exp(+II*k);
	Hk[fq(1,0,8)] += -tIO*exp(+II*(-Ax_t(time)-sqrt(3.)*Ay_t(time))/2.)*exp(+II*k/2.);
	Hk[fq(7,1,8)] += -tIO*exp(+II*(-Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(+II*k/2.);
	Hk[fq(6,1,8)] += -tI2*exp(-II*Ax_t(time))*exp(+II*k);
	Hk[fq(2,1,8)] += -tI1*exp(+II*(-Ax_t(time)-sqrt(3.)*Ay_t(time))/2.)*exp(+II*k/2.);
	Hk[fq(6,2,8)] += -tI1d*exp(+II*(-Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(+II*k/2.);
	Hk[fq(5,2,8)] += -tI2*exp(-II*Ax_t(time))*exp(+II*k);
	Hk[fq(3,2,8)] += -tIO*exp(+II*(-Ax_t(time)-sqrt(3.)*Ay_t(time))/2.)*exp(+II*k/2.);
	Hk[fq(5,3,8)] += -tIO*exp(+II*(-Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(+II*k/2.);
	Hk[fq(4,3,8)] += -tOd*exp(-II*Ax_t(time))*exp(+II*k);
	Hk[fq(5,4,8)] += -tIO*exp(+II*(+Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(-II*k/2.);
	Hk[fq(6,5,8)] += -tI1*exp(+II*(+Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(-II*k/2.);
	Hk[fq(7,6,8)] += -tIO*exp(+II*(+Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(-II*k/2.);
	Hk[fq(4,3,8)] += -tO*exp(+II*Ax_t(time))*exp(-II*1.*k);
	Hk[fq(4,2,8)] += -tIO*exp(-II*(-Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(-II*k/2.);
	Hk[fq(5,2,8)] += -tI2*exp(+II*Ax_t(time))*exp(-II*1.*k);
	Hk[fq(5,1,8)] += -tI1d*exp(-II*(-Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(-II*k/2.);	
	Hk[fq(6,1,8)] += -tI2*exp(+II*Ax_t(time))*exp(-II*1.*k);	
	Hk[fq(6,0,8)] += -tIO*exp(-II*(-Ax_t(time)+sqrt(3.)*Ay_t(time))/2.)*exp(-II*k/2.);	
	Hk[fq(7,0,8)] += -tO*exp(+II*Ax_t(time))*exp(-II*1.*k);	
	
	// Upper triagonal
	Hk[fq(0,7,8)] = conj(Hk[fq(7,0,8)]);
	Hk[fq(0,1,8)] = conj(Hk[fq(1,0,8)]);
	Hk[fq(1,7,8)] = conj(Hk[fq(7,1,8)]);
	Hk[fq(1,6,8)] = conj(Hk[fq(6,1,8)]);
	Hk[fq(1,2,8)] = conj(Hk[fq(2,1,8)]);
	Hk[fq(2,6,8)] = conj(Hk[fq(6,2,8)]);
	Hk[fq(2,5,8)] = conj(Hk[fq(5,2,8)]);
	Hk[fq(2,3,8)] = conj(Hk[fq(3,2,8)]);
	Hk[fq(3,5,8)] = conj(Hk[fq(5,3,8)]);
	Hk[fq(3,4,8)] = conj(Hk[fq(4,3,8)]);
	Hk[fq(4,5,8)] = conj(Hk[fq(5,4,8)]);
	Hk[fq(5,6,8)] = conj(Hk[fq(6,5,8)]);
	Hk[fq(6,7,8)] = conj(Hk[fq(7,6,8)]);
	Hk[fq(3,4,8)] = conj(Hk[fq(4,3,8)]);
	Hk[fq(2,4,8)] = conj(Hk[fq(4,2,8)]);
	Hk[fq(2,5,8)] = conj(Hk[fq(5,2,8)]);
	Hk[fq(1,5,8)] = conj(Hk[fq(5,1,8)]);
	Hk[fq(1,6,8)] = conj(Hk[fq(6,1,8)]);
	Hk[fq(0,6,8)] = conj(Hk[fq(6,0,8)]);
	Hk[fq(0,7,8)] = conj(Hk[fq(7,0,8)]);
}


void groundstate(cvec &Hk, dvec &evals, dvec &BZ, double &mu, int &numprocs, int &myrank)
/**
 *	Calculation of chemical potential
 *  -Hk: Complex vector[64] to store Hamiltonian
 *  -evals: Real vector[8] of eigenvalues
 *  -BZ: k-points of reduced 1d reciprocal cell
 *  -mu: Chemical potential
 *  -numprocs: Total number of processes (MPI)
 *  -myrank: Rank of process (MPI)
 */
{	
	int count = 0;                                                      // count # of loops of self-consistency
	double N_tot;	
	double mu_old;    
	double deviation = 1.0;
	mu = mu_init;
	
	vector<dvec> N0(NN, dvec(8));                                       // Density in band basis	
	while(deviation > dev)
	{
		count++;
					
		mu_old = mu;	
	    
		N_tot = 0.;
		for(int k=myrank; k<NN; k+=numprocs)
		{		
			set_Hk(BZ[k], Hk, 0.0);
			diagonalize(Hk, evals);                                 	// Hk -> S (eigenvectors as columns)
			for(int i=0; i<8; i++)
			{			
				N0[k][i] = fermi(evals[i], mu);
				N_tot +=  N0[k][i];	
			}
		}

#ifndef NO_MPI		
		MPI_Allreduce(MPI_IN_PLACE, &N_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif		
		mu += -DELTA*(N_tot-4.*double(NN));	                       // constant value method to get new mu
		
		deviation = abs(mu-mu_old);
		if(myrank==0){
			cout << "loop #" << count << ": deviation = " << deviation << endl;
			cout << "chemical potential mu = " << mu << endl;
			cout << "Particle number per unit cell = " << N_tot/double(NN) << endl; 
		}	
	}
	if(myrank==0)
	{
		ofstream myfile ("mu.txt");
		if (myfile.is_open())
		{
			myfile << mu;
			myfile.close();
		}	
		else cout << "Unable to open file" << endl;	
	}
}

	
void Hk_bands(vector<dvec> &BANDS, cvec &Hk, dvec &evals, dvec &BZ, const string& filename)
/**
 *	Calculate bands of Hk(k) for high symmetry path and store them in "bands.txt":
 *	-BANDS: Vector of real vectors[3] to store k-points of path
 *  -Hk: Complex vector[64] to store Hamiltonian
 *  -evals: Real vector[8] of eigenvalues
 *  -BZ: k-points of reduced 1d reciprocal cell
 *	-filename: String to define file
 */
{
	for(int k=0; k<NN; k++)
	{
		set_Hk(BZ[k], Hk, 0.0);
		diagonalize(Hk, evals);
		for(int m=0; m<8; m++)
			BANDS[k][m] = evals[m];
	}
	ofstream myfile (filename);
	if (myfile.is_open())
	{
		for(int k=0; k<NN; k++)
		{
			for(int m=0; m<8; m++)
			{
				myfile << BANDS[k][m] << " " ;
			}
		myfile  << endl;
		}
	myfile.close();
	}
    else cout << "Unable to open file" << endl;
}


void Propagation(cvec &Hk, dvec &evals, vector<cvec> &UMATRIX, double &mu, dvec &BZ, vector<cvec> &RHO_t, dvec ETOT_t, int &numprocs, int &myrank)
/**
 *	Mid-point Euler propagation scheme for density matrix:
 *  -Hk: Complex vector[64] to store Hamiltonian
 *  -evals: Real vector[8] of eigenvalues
 *  -UMATRIX: Vector[TIMESTPES] of complex vectors[64] to store propagators 
 * 	-mu: chemical potential
 *  -BZ: k-points of reduced 1d reciprocal cell
 *	-RHO_t: Vector[TIMESTPES] of complex vector[64] to store t.-d. density matrix
 *  -ETOT_t: Real vector[TIMESTEPS] to store total t.-d. total energy
 *	-numprocs: Total number of processes (MPI)
 *	-myrank: Rank of process (MPI) 
 */
{
	double h = (endtime-starttime)/TIMESTEPS;
		
	cvec TEMP1(64); 
	cvec TEMP2(64); 
			
	for(int t=0; t<TIMESTEPS; t++)
		ETOT_t[t] = 0.0;
	
	// Calculation of mid-point Euler propagators U(t,t+dt)
	for(int k=myrank; k<NN; k+=numprocs)
	{
		if(myrank==0) cout << "k = " << k << endl;
		for(int t=0; t<TIMESTEPS-1; t++)
		{
			set_Hk(BZ[k], TEMP1, h*double(t));	
			set_Hk(BZ[k], TEMP2, h*double(t+1));
			for(int i=0; i<64; i++)
				Hk[i] = 0.5*(TEMP1[i]+TEMP2[i]);
			diagonalize(Hk, evals); 
			for(int i=0; i<8; i++)
			{
				for(int j=0; j<8; j++)
				{
					TEMP1[fq(i,j,8)] = exp(+II*evals[i]*h)*delta(i,j);
				}	
			}
			times(TEMP1, Hk, TEMP2);                                        // S @ H @ S^(-1) = H_D --> nk = S^(-1) @ nk_D @ S
			times_dn(Hk, TEMP2, UMATRIX[t]);		
		}	
		// Set initial Density in band basis
		set_Hk(BZ[k], Hk, starttime);
		diagonalize(Hk, evals);
		for(int i=0; i<8; i++)
		{
			for(int j=0; j<8; j++)
			{
				RHO_t[0][fq(i,j,8)] = fermi(evals[i], mu)*delta(i,j);
			}	
		}
		times(RHO_t[0], Hk, TEMP1);                                        // S @ H @ S^(-1) = H_D --> nk = S^(-1) @ nk_D @ S
		times_dn(Hk, TEMP1, RHO_t[0]);	

		// Popagation of initial density operator and calculation of total energy
		for(int t=0; t<TIMESTEPS-1; t++)
		{	
			times(RHO_t[t], UMATRIX[t], TEMP1);
			times_dn(UMATRIX[t], TEMP1, RHO_t[t+1]);
			set_Hk(BZ[k], Hk, h*double(t));	                             
			times(RHO_t[t], Hk, TEMP1);
			for(int i=0; i<8; i++)
			{
				ETOT_t[t] += real(TEMP1[fq(i,i,8)])/double(NN);	
			}			
		}
		ETOT_t[TIMESTEPS-1] = ETOT_t[TIMESTEPS-2];
	}
#ifndef NO_MPI		
	MPI_Allreduce(MPI_IN_PLACE, &ETOT_t[0], TIMESTEPS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif	
	if(myrank==0)
	{
		const string filename = "ETOT_t.txt";
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			for(int t=0; t<TIMESTEPS; t++)
			{
				// Add factor 2 for Spin degeneracy
				myfile << 2.0*ETOT_t[t] << endl;                        
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
}	 


// main() function #####################################################

int main(int argc, char * argv[])
{
    //************** MPI INIT ***************************
  	int numprocs=1, myrank=0, namelen;
    
#ifndef NO_MPI
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Get_processor_name(processor_name, &namelen);
    
	cout << "Process " << myrank << " on " << processor_name << " out of " << numprocs << " says hello." << endl;
	MPI_Barrier(MPI_COMM_WORLD);
    
#endif
	if(myrank==0) cout << "\n\tProgram running on " << numprocs << " processors." << endl;

	// DECLARATION AND INTITALIZATIO
	// chemical potential
	double mu = 0.0;
	
	//Sampling of reduiced Brillouin zone [-pi/2,+pi/2] 
	dvec BZ(NN);
	for(int ii=0; ii<NN; ii++)
	{
		BZ[ii] = -PI/2. + PI*double(ii)/double(NN);                        
	}
	
	// allocation of matrices RHO[k]
	vector<cvec> RHO_0(NN, cvec(8,0.0));  
	
	// vector for eigenvalues
	dvec evals(8);

	// vector for Hamiltonian Hk
	cvec Hk(64);
	
	// bands 
	vector<dvec> BANDS(NN,dvec(8));
	
	// Vector to store t.-d. total energy
	dvec ETOT_t(TIMESTEPS);
	
	// Store Euler propagators
	vector<cvec> UMATRIX(TIMESTEPS,cvec(64));      
	
	// Vector to store t.-d. density operator
	vector<cvec> RHO_t(TIMESTEPS, cvec(64));                                              
	
	// CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	const clock_t begin_time = clock();	
	
	double h = (endtime-starttime)/TIMESTEPS;
	
	// Store t.-d. gauge field in file
	double Ax;
	double Ay;
	if(myrank==0)
	{
		ofstream myfile ("DRIVING_t.txt");
		if (myfile.is_open())
		{
			myfile << "#Nuclear"  << " " << "Peierls" << " " << endl;
			for(int t=0; t<TIMESTEPS; t++)
			{			
				Ax = Ax_t(h*double(t));
				Ay = Ay_t(h*double(t));
				myfile << Ax << " " << Ay << " " << endl;	
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	
	// Calculate initial density and mu
	groundstate(Hk, evals, BZ, mu, numprocs, myrank);
	
	// Calculate equilibrium band structure
	if(myrank==0) 
	{
		//calulation of initial bands
		Hk_bands(BANDS, Hk, evals, BZ, "bands0.txt");
	}
	
	// Popagate and get toal energy
	cout << Ax_peierls << "#######################################################################################" << endl; 
	Propagation(Hk, evals, UMATRIX, mu, BZ, RHO_t, ETOT_t, numprocs, myrank);
	
	if(myrank==0) cout << "Calculations lasted: " << float(clock() - begin_time)/CLOCKS_PER_SEC << " seconds" << endl;
	
#ifndef NO_MPI
	MPI_Finalize();
#endif		
}

