/***********************************************************
* Copyright (C) 2014  
* Authors: Hamid Teimouri & Daniel Celis
* Rice university--Department of Chemistry
* This file is distributed under the terms of the
* GNU General Public License as published by the
* Free Software Foundation; either version 3 of the
* License, or (at your option) any later version.
* http://www.gnu.org/copyleft/gpl.txt
***********************************************************/
//===========================================================================================================//
// Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP) with interacting particles //
//================================= for periodic boundary conditions (ring) =================================//
//===========================================================================================================//

// For details on the exact solution see:
//J. S. Hager, J. Krug, V. Popkov, and G. M. Schuts, "Minimal current phase and universal boundary layers in driven diffusive systems." Physical Review E, Volume 63, 056110.

#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <math.h>
#include <time.h>
#include <numeric>
#include <cstdlib>
#include <vector>
#include <valarray>
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <ctime>
#include <cmath>
#include "ran3.h"
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

using std::vector;
using namespace std;

long int dum;

const int L = 100;
long double T = 1e4;
long double Teq = T * 0.20;
long double Tdif = T - Teq;
long double dt = 0.01;
long double t;

int site;
int prevsite;
int nextsite;
int nnextsite;

long double J;
long double SiteDens;
long double AvgDens;

long double rho;
long double epsilon;
long double theta;
long double q;
long double r;
long double densprof[L+1], lattice[L+1], p[L+1], possibleEnter[L+1];
long double n_i_hop, q_hop, r_hop;

long int j, num, k, i, N;

double lambda;
double Jtheory;

///////////////////////// Functions.////////////////////////
void possible_enter()
{
	i = 1;
	k = 1;
	while (i < L)
	{
		if (lattice[i] == 0)
		{
			possibleEnter[k] = i;
			i++;
			k++;
		}
		else
		{
			i++;
		}	
	}
}

void initialise() // Initialize lattice.
{
	t = 0.0; // Initial time.
	J = 0.0; // Initial current
	AvgDens = 0.0;
	num = 0;
	for (j = 1; j <= L; j++)
	{
		lattice[j]=0; densprof[j]=0;  p[j]=1; possibleEnter[j]=0;
	}
	n_i_hop = dt * p[L/2]; // The probability of hopping.
	q_hop = q * n_i_hop; // The probability of hopping.
	r_hop = r * n_i_hop; // The probability of hopping.
	while (num < N)
	{	
		possible_enter();
		site = possibleEnter[rand()%k + 1];
		lattice[site] = 1;
		num++;
	}
}

void boundary_interactions()
{ 	
	if (site == 1 && lattice[site] == 1 && lattice[nextsite] == 0)
	{
		// Ejection rate beta if in state (0,1)
		if (lattice[L] == lattice[nnextsite] && ran3(&dum) <= n_i_hop)
		{
			lattice[site] = 0;
			lattice[nextsite] = 1;
		}
		// Ejection rate r*beta if in state (1,1) => Disassociation
		if (lattice[L] != lattice[nnextsite])
		{
			if(lattice[L] == 1 && ran3(&dum) <= r_hop)
			{
				lattice[site] = 0;
				lattice[nextsite] = 1;
			}
			if(lattice[L] == 0 && ran3(&dum) <= q_hop)
			{
				lattice[site] = 0;
				lattice[nextsite] = 1;
			}
			
		}
	}
	
	if (site == L && lattice[site] == 1 && lattice[1] == 0)
	{
		// Ejection rate beta if in state (0,1)
		if (lattice[L-1] == lattice[2] && ran3(&dum) <= n_i_hop)
		{
			lattice[site] = 0;
			lattice[1] = 1;
		}
		// Ejection rate r*beta if in state (1,1) => Disassociation
		if (lattice[L-1] != lattice[2])
		{
			if(lattice[2] == 1 && ran3(&dum) <= q_hop)
			{
				lattice[site] = 0;
				lattice[1] = 1;
			}
			if(lattice[2] == 0 && ran3(&dum) <= r_hop)
			{
				lattice[site] = 0;
				lattice[1] = 1;
			}
			
		}
	}
	if (site == L-1 && lattice[site] == 1 && lattice[L] == 0)
	{
		// Ejection rate beta if in state (0,1)
		if (lattice[L-2] == lattice[1] && ran3(&dum) <= n_i_hop)
		{
			lattice[site] = 0;
			lattice[L] = 1;
		}
		// Ejection rate r*beta if in state (1,1) => Disassociation
		if (lattice[L-2] != lattice[1])
		{
			if(lattice[1] == 1 && ran3(&dum) <= q_hop)
			{
				lattice[site] = 0;
				lattice[L] = 1;
			}
			if(lattice[1] == 0 && ran3(&dum) <= r_hop)
			{
				lattice[site] = 0;
				lattice[L] = 1;
			}
			
		}
	}		
}

void move()
{
	for (j = 1; j <= L; j++)
	{
		site = rand()%L + 1; 
		prevsite = site - 1;
		nextsite = site + 1;
		nnextsite = nextsite + 1;
		boundary_interactions();
		if (lattice[site] == 1 && lattice[nextsite] == 0 && site >= 2 && site <= L-2)
		{
			if (lattice[prevsite] == lattice[nnextsite] && ran3(&dum) <= n_i_hop) // Bulk hopping with rate 1 if in state (0,1,0,0) or (1,1,0,1)
			{ 
				lattice[nextsite] = 1;      // moves to state (0,0,1,0) 
				lattice[site] = 0; 
				if (site == L/2 && t >= Teq)
				{
					
					J++;	
				}
			}
			
			if(lattice[prevsite] != lattice[nnextsite])
			{
				if(lattice[prevsite] == 1 && ran3(&dum) <= r_hop) // Bulk hopping with rate 1 if in state (1,1,0,0)
				{  
					lattice[nextsite] = 1;      // moves to state (1,0,1,0)
					lattice[site] = 0;      
					if (site == L/2 && t >= Teq)
					{
						J++;
					}
				}
				if(lattice[prevsite] == 0 && ran3(&dum) <= q_hop) // Bulk hopping with rate q if in state (0,1,0,1) => Binding
				{
					lattice[nextsite] = 1;      // moves to state (0,0,1,1) 
					lattice[site] = 0;
					if (site == L/2 && t >= Teq)
					{
						J++;
					}	
				}
			}
		} // End of site availability check.
	} // End of loop through array.
}

void update()
{
	dum=-time(NULL);
	ran3(&dum);
	for (t = 0.0; t <= T; t += dt) // Time loop allows for better averaging, and allows the system to reach steady state.
	{
		move(); // Move function.
		if (t >= Teq) // Building the density profile after the system has reached steady state. Not averaged yet.
		{
			for (j = 1; j <= L; j++) // Sweeps throught the array.
			{
				if (lattice[j] == 1) // If it finds a site which contains a particle.
				{
					densprof[j] += dt; // It adds one to its density counter. In reality this can be >1, it will be time-averaged later.
				}
			}
		}
	} // End of time loop.
}

double cputime ( )
{

  double time;
  time = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;
  
  return time;
}
//===================================================================//
//============================ MAIN CODE ============================//
//===================================================================//
int main()
{
	double cputime0;
	double cputime1;
	double cputime2;  
	const string program ="Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP)";
	const string spaces(program.size(), '*');
	const string stars = spaces;
	cout<<"\n"<<endl;
	cout<<stars<<endl;
	cout<< program <<endl;
	cout << "for interacting particles with periodic boundaries (ring)." <<endl;
	cout<<stars<<endl;
	
	char key[] = "n";
	char buffer[80];
	string tc;
	
	cout<< "\nrho = ";
	cin>> rho;
	cin.ignore();
	
	N = int(L*rho);
	
	printf ("Thermodynamically consistent attraction and repulsion rates (y/n)? Default: yes. If no, 0 <= epsilon <= 1.\n");
	fflush (stdout);
	scanf ("%79s",buffer);
	
	cout<< "E = ";
	cin>> epsilon;
	cin.ignore();
	
	if (strcmp (key,buffer) != 0)
	{
		cout<< "theta = ";
		cin>> theta;
		cin.ignore();
		
		q=exp(theta*epsilon);
		r=exp((theta-1)*epsilon);
		
		tc = "C";
		
		if (epsilon == 0.0)
		{
			cout<< "\n" << stars <<endl;
			cout<< "Exact Solution:\n" <<endl;
			cout<< " J = rho(1-rho) = "<< rho*(1-rho) <<endl;
			cout<< " Density = rho = "<< rho <<endl;
			cout<< "\n" << stars <<endl;
		}		
	}
	
	if (strcmp (key,buffer) == 0)
	{
		q = 1 - epsilon;
		r = 1 + epsilon;
		
		lambda = 1 / ( sqrt(4 * rho * ( 1 - rho ) ) ) + sqrt( 1/ (4 * rho * ( 1 - rho ) ) - 1 + q/r );
		Jtheory = ( lambda - epsilon * sqrt( 4 * rho * ( 1 - rho ) ) )/ pow( lambda, 3 );
		
		cout<< "\n" << stars <<endl;
		cout<< "Exact Solution:\n" <<endl;
		cout<< " J. S. Hager, J. Krug, V. Popkov, and G. M. Schuts" << endl;
		cout<< " ''Minimal current phase and universal boundary layers in driven diffusive systems.''" <<endl;
		cout<< " Physical Review E, Volume 63, 056110.\n" <<endl;
		cout<< " lambda = 1/( sqrt( 4*rho*( 1-rho ) ) ) + sqrt( 1/( 4*rho*( 1-rho ) ) - 1 + q/r ) = " << lambda <<endl;
		cout<< " J = ( lambda - epsilon * sqrt( 4*rho*( 1-rho ) ) )/pow( lambda, 3 ) = " << Jtheory <<endl;
		cout<< " Density = rho = "<< rho <<endl;
		cout<< "\n" << stars <<endl;
		
		tc = "I";
	}
	cout<< "Parameters:\n" <<endl;
	cout<< " L = " << L << endl;
	cout<< " T = " << T << endl;
	cout<< " dt = " << dt <<endl;
	cout<< " rho = " << rho <<endl;
	cout<< " N = " << N <<endl;
	cout<< " E = " << epsilon << endl;
	if (strcmp (key,buffer) != 0)
	{
		cout<< " theta = " << theta << endl;
	}
	cout<< " q = " << q << endl;
	cout<< " r = " << r << endl;
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamD("");
	if (strcmp (key,buffer) != 0)
	{
		fileNameStreamD << tc << "_R_D" << "_E=" << epsilon << "_th=" << theta <<"_rho="<< rho <<".txt";
	}
	if (strcmp (key,buffer) == 0)
	{
		fileNameStreamD << tc << "_R_D" << "_E=" << epsilon <<"_rho="<< rho <<".txt";
	}
    	std::string fileNameD = fileNameStreamD.str();
	ofstream q1(fileNameD.c_str());
	
	#pragma omp parallel
	#pragma omp for lastprivate(lattice, J)
	initialise();
	update();
	
	for (j = 1; j <= L; j++) // Time averaging of the density.
	{
		SiteDens = densprof[j]/Tdif;
		q1<<j<<" "<<SiteDens<<endl;
		AvgDens += SiteDens;
	}
  
	q1.close();
        cputime2 = cputime ();
        cputime0 = cputime2 - cputime1;
	cout<< "Simulation Results:\n" <<endl;
	cout<<" J = "<< J/Tdif <<endl;
	cout<<" Bulk Density = "<< AvgDens/L <<endl;
	cout<< "\n";
	cout<< " Elapsed cpu time for main computation:\n";
	cout<< "  " << cputime2 << " seconds" <<endl;
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamL("");
	if (strcmp (key,buffer) != 0)
	{
    		fileNameStreamL << tc << "_R_L" << "_E=" << epsilon << "_th=" << theta <<"_rho="<< rho <<".txt";
	}
	if (strcmp (key,buffer) == 0)
	{
		fileNameStreamL << tc << "_R_L" << "_E=" << epsilon <<"_rho="<< rho <<".txt";
	}
	std::string fileNameL = fileNameStreamL.str();
	ofstream q2(fileNameL.c_str());
	
	if (strcmp (key,buffer) != 0)
	{
		q2<< "\tThermodynamically consistent.\n" <<endl;
	}
	if (strcmp (key,buffer) == 0)
	{
		q2<< "\tThermodynamically inconsistent.\n" <<endl;
	}
	q2<< stars <<endl;
	q2<< "Parameters:\n" <<endl;
	q2<< " L = " << L << endl;
	q2<< " T = " << T << endl;
	q2<< " dt = " << dt <<endl;
	q2<< " rho = " << rho <<endl;
	q2<< " N = " << N <<endl;
	q2<< " E = " << epsilon << endl;
	if (strcmp (key,buffer) != 0)
	{
		q2<< " theta = " << theta << endl;
	}
	q2<< " q = " << q << endl;
	q2<< " r = " << r << endl;
	q2<< "\n" << stars <<endl;
	
	if (strcmp (key,buffer) != 0)
	{
		
		if (epsilon == 0.0)
		{
			q2<< "Exact Solution:\n" <<endl;
			q2<< " J = rho(1-rho) = "<< rho*(1-rho) <<endl;
			q2<< " Density = rho = "<< rho <<endl;
		}		
	}
	
	if (strcmp (key,buffer) == 0)
	{	
		q2<< "Exact Solution:\n" <<endl;
		q2<< " J. S. Hager, J. Krug, V. Popkov, and G. M. Schuts" << endl;
		q2<< " ''Minimal current phase and universal boundary layers in driven diffusive systems.''" <<endl;
		q2<< " Physical Review E, Volume 63, 056110.\n" <<endl;
		q2<< " lambda = 1/( sqrt( 4*rho*( 1-rho ) ) ) + sqrt( 1/( 4*rho*( 1-rho ) ) - 1 + q/r ) = " << lambda <<endl;
		q2<< " J = ( lambda - epsilon * sqrt( 4*rho*( 1-rho ) ) )/pow( lambda, 3 ) = " << Jtheory <<endl;
		q2<< " Density = rho = "<< rho <<endl;
		q2<< "\n" << stars <<endl;
	}
	
	q2<< "Simulation Results:\n" <<endl;
	q2<<" J = "<< J/Tdif <<endl;
	q2<<" Bulk Density = "<<AvgDens/L<<endl;
	q2<< "\n";
	q2<< " Elapsed cpu time for main computation:\n";
	q2<< "  " << cputime2 << " seconds" <<endl;
	q2<< "\n" << stars <<endl;
	q2.close();
}
