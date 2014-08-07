/**********************************************************************************************************
* Copyright (C) 2014  
* Authors: Hamid Teimouri & Daniel Celis
* Rice university--Department of Chemistry
* This file is distributed under the terms of the
* GNU General Public License as published by the
* Free Software Foundation; either version 3 of the
* License, or (at your option) any later version.
* http://www.gnu.org/copyleft/gpl.txt
************************************************************************************************************/
//===========================================================================================================//
// Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP) with interacting particles //
//===========================================================================================================//

// For details on the exacto solution see:
// B. Derrida, M.R. Evans, V. Hakim, V. Pasquier, Exact solution of a 1d asymmetric exclusion model using a matrix formulation J. Phys. A26, 1493-1517 (1993)

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
long double Teq = T * 0.2;
long double Tdif = T - Teq;
long double dt = 0.1;
long double t;

int site;
int prevsite;
int nextsite;
int nnextsite;

long double J;
long double SiteDens;
long double AvgDens;

long double alpha;
long double beta;
long double epsilon;
long double theta;
long double q;
long double r;
long double densprof[L+1], lattice[L+2], p[L+1];
long double n_i_enter, i_enter, n_i_eject, i_eject, n_i_hop, q_hop, r_hop;

int j;

///////////////////////// Functions.////////////////////////
void initialise()
{
	J = 0.0;
	AvgDens = 0.0;
	lattice[0] = 0;
	lattice[L+1] = 0;
	for (j = 1; j <= L; j++)
	{
		lattice[j]=0; densprof[j]=0;  p[j]=1;
	}
	n_i_enter = dt*alpha*p[1];
	i_enter = q*n_i_enter;
	n_i_eject = dt*beta*p[L];
	i_eject = r*n_i_eject;
	n_i_hop = dt*p[L/2];
	q_hop = dt*q*p[L/2];
	r_hop = dt*r*p[L/2];	
}
void boundary_interactions()
{ 
	if (site == 1 && lattice[site] == 0)
	{
		// Injection with rate alpha if in state (0,0)
		if (lattice[nextsite] == 0 && ran3(&dum) <= n_i_enter)
		{
			lattice[site] = 1;
		}
		// Injection with rate q*alpha if in state (0,1) => Binding
		if (lattice[nextsite] == 1 && ran3(&dum) <= i_enter)
		{
			lattice[site] = 1;
		}
	}	
	if (site == L && lattice[site] == 1)
	{
		// Ejection rate beta if in state (0,1)
		if (lattice[prevsite] == 0 && ran3(&dum) <= n_i_eject)
		{
			lattice[site] = 0;
		}
		// Ejection rate r*beta if in state (1,1) => Disassociation
		if (lattice[prevsite] == 1 && ran3(&dum) <= i_eject)
		{
			lattice[site] = 0;
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
		if (lattice[site] == 1 && lattice[nextsite] == 0 && site >= 1 && site <= L-1)
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
			if (lattice[prevsite] != lattice[nnextsite])
			{
				if (lattice[prevsite] == 1 && ran3(&dum) <= r_hop)
				{
					lattice[nextsite] = 1;      // moves to state (1,0,1,0)
					lattice[site] = 0;      
					if (site == L/2 && t >= Teq)
					{
						J++;
					}
				}
				if (lattice[prevsite] == 0 && ran3(&dum) <= q_hop)
				{
					lattice[nextsite] = 1;      // moves to state (1,0,1,0)
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
	const string program ="Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP) with interacting particles.";
	const string spaces(program.size(), '*');
	const string stars = spaces;
	cout<<"\n"<<endl;
	cout<<stars<<endl;
	cout<< program <<endl;
	cout<<stars<<endl;
	
       //ask user for entrace and exit rates
	cout<< "\nalpha = ";
	cin>> alpha;
	cin.ignore();
	cout<< "beta = ";
	cin>> beta;
	cin.ignore();
         
       //ask user for interaction energy
	cout<< "E = ";
	cin>> epsilon;
	cin.ignore();
	//ask user for the interaction energy spliting fraction
	if (epsilon != 0.0)
	{
		cout<< "theta = ";
		cin>> theta;
		cin.ignore();
	}

        // binding & unbinding rates are related in a thermodynamically consistent fashion
	q = exp(theta * epsilon); // binding rate
	r = exp((theta - 1) * epsilon); // unbinding rate
        
	cout<< "\n" << stars <<endl;
	cout<< "Parameters:\n" <<endl;
	cout<< " L = " << L << endl;
	cout<< " T = " << T << endl;
	cout<< " dt = " << dt <<endl;
	cout<< " alpha = " << alpha << endl;
	cout<< " beta = " << beta << endl;
	cout<< " E = " << epsilon << endl;
	cout<< " theta = " << theta << endl;
	cout<< " q = " << q << endl;
	cout<< " r = " << r << endl;
	cout<< "\n" << stars <<endl;

	if (epsilon == 0.0)
	{
		cout<< "Exact Solution:\n" <<endl;
		cout<< "B. Derrida, M.R. Evans, V. Hakim, V. Pasquier, Exact solution of a 1d asymmetric exclusion model using a matrix formulation J. Phys. A26, 1493-1517 (1993)" <<endl;
		if (alpha > 0.5 && beta > 0.5)
		{
			cout<< " Maximal Current Phase (q=r=1)" <<endl;
			cout<< "\n";
			cout<< " J = 0.25" <<endl;
			cout<< " Bulk Density = 1/2" <<endl;
		}	
		if (beta < alpha && beta < 0.5)
		{	
			cout<< " High Density Phase (q=r=1)" <<endl;
			cout<< "\n";
			cout<< " J = b(1-b) = " << beta*(1-beta) <<endl;
			cout<< " Bulk Density = " << 1.0-beta <<endl;
		}	
		if (alpha < beta && alpha < 0.5)
		{	
			cout<< " Low Density Phase (q=r=1)"<<endl;
			cout<< "\n";
			cout<< " J = a(1-a) = " << alpha*(1-alpha) <<endl;
			cout<< " Bulk Density = " << alpha <<endl;
		}
		cout<< "\n" << stars <<endl;
	}
	
	std::ostringstream fileNameStreamD("");
    	fileNameStreamD << "Densprof" << "_E=" << epsilon << "_th=" << theta <<"_a=" << alpha << "_b="<< beta <<".txt";
    	std::string fileNameD = fileNameStreamD.str();
    	ofstream q1(fileNameD.c_str());
	
	#pragma omp parallel
	#pragma omp for lastprivate(lattice, J)  
	
	initialise();
	update();
	
	for (j = 1;j <= L; j++) // Time averaging of the density. It sweeps through the array.
	{
		SiteDens = densprof[j]/Tdif;
		q1<< j <<" "<< SiteDens <<endl;
		AvgDens += SiteDens;
	}
	q1.close();
	cputime2 = cputime ();
	cputime0 = cputime2 - cputime1;
	cout<< "Simulation Results:\n" <<endl;
	cout<<" J = "<< J/Tdif <<endl;
	cout<< " Bulk Density = " << AvgDens/L <<endl;
	cout<< "\n";
	cout<< " Elapsed cpu time for main computation:\n";
	cout<< "  " << cputime2 << " seconds";
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamL("");
    	fileNameStreamL << "Log" << "_E=" << epsilon << "_th=" << theta <<"_a=" << alpha << "_b="<< beta <<".txt";
    	std::string fileNameL = fileNameStreamL.str();
    	ofstream q2(fileNameL.c_str());
	
	q2<< stars <<endl;
	q2<< "Parameters:\n" <<endl;
	q2<< " L = " << L << endl;
	q2<< " T = " << T << endl;
	q2<< " dt = " << dt <<endl;
	q2<< " alpha = " << alpha << endl;
	q2<< " beta = " << beta << endl;
	q2<< " E = " << epsilon << endl;
	q2<< " theta = " << theta << endl;
	q2<< " q = " << q << endl;
	q2<< " r = " << r << endl;
	q2<< "\n" << stars << endl;
	
	if (epsilon == 0.0)
	{
		q2<< "Exact Solution:\n" <<endl;
		q2<< "B. Derrida, M.R. Evans, V. Hakim, V. Pasquier, Exact solution of a 1d asymmetric exclusion model using a matrix formulation J. Phys. A26, 1493-1517 (1993)" <<endl;
		if (alpha > 0.5 && beta > 0.5)
		{
			q2<< " Maximal Current Phase (q=r=1)" <<endl;
			q2<< "\n";
			q2<< " J = 0.25" <<endl;
			q2<< " Bulk Density = 1/2" <<endl;
		}	
		if (beta < alpha && beta < 0.5)
		{	
			q2<< " High Density Phase (q=r=1)" <<endl;
			q2<< "\n";
			q2<< " J = b(1-b) = " << beta*(1-beta) <<endl;
			q2<< " Bulk Density = " << 1.0-beta <<endl;
		}	
		if (alpha < beta && alpha < 0.5)
		{	
			q2<< " Low Density Phase (q=r=1)"<<endl;
			q2<< "\n";
			q2<< " J = a(1-a) = " << alpha*(1-alpha) <<endl;
			q2<< " Bulk Density = " << alpha <<endl;
		}
		q2<< "\n" << stars <<endl;
	}
	
	q2<< "Simulation Results:\n" <<endl;
	q2<< " J = "<< J/Tdif <<endl;
	q2<< " Bulk Density = "<< AvgDens/L <<endl;
	q2<< "\n";
	q2<< " Elapsed cpu time for main computation:\n";
	q2<< "  " << cputime2 << " seconds" <<endl;
	q2<< "\n" << stars <<endl;
	q2.close();
}
