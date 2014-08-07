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
//==============================================================================================================//
// Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP) for periodic boundaries (ring)//
//==============================================================================================================//

#include <iostream>
#include <fstream>
#include <sstream>
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
#include "ran3.h"
#include <cmath>
#pragma hdrstop
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
long double Tdif = T-Teq;
long double dt = 0.01; 
long double t;

int site;
int nextsite;

long double J;
long double SiteDens;
long double AvgDens;


long double hop;
long double rho;
double densprof[L+1], lattice[L+1], p[L+1], possibleEnter[L+1];

long int j, num, k, i, N;

//////////////////////////////////////// Functions ///////////////////////////////////////
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
	t=0.0; // Initial time.
	J=0.0; // Initial current
	AvgDens = 0.0;
	num = 0;
	for (j = 1; j <= L; j++)
	{
		lattice[j]=0; densprof[j]=0;  p[j]=1;
	}
	hop = dt * p[L/2]; // The probability of hopping.
	while (num < N)
	{	
		possible_enter();
		site = possibleEnter[rand()%k + 1];
		lattice[site] = 1;
		num++;
	}
}

void move()
{	
	for (j = 1;j <= L; j++) // Goes through as many iterations as the array's length.
	{
		site = rand()%L + 1; // Picks out a random site in the array.
		nextsite = site + 1; // Defines the next site. Saves the program from unecessary calculations.
		if (site == L && lattice[site] == 1 && lattice[1] == 0 && ran3(&dum) <= hop) // Exit last site.
		{
			lattice[site] = 0; // If the citeria are met, the site is emptied.
			lattice[1] = 1;
		}
		if (lattice[site] == 1 && lattice[nextsite] == 0 && site < L && ran3(&dum) <= hop) // Bulk of the array.
		{
			lattice[site] = 0; // If the criteria are met the particle leave site.
			lattice[nextsite] = 1; // And enters the next site.
			if(site == L/2 && t >= Teq) // Measuring current. It's measured when a particle crosses the midpoint of the array.
			{
				J++;
			}	
		}		
	
	} // End of loop through array. J-loop.
	
}

void update()
{
	dum=-time(NULL);
	ran3(&dum);
	for (t = 0.0; t <= T; t += dt)
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
	cout<< "for periodic boundary conditions (ring)." <<endl;
	cout<<stars<<endl;
	
	cout<< "\nrho = ";
	cin>> rho;
	cin.ignore();
	
	N = int(L*rho);
	
	cout<< "\n" << stars <<endl;
	cout<< "Parameters:\n" <<endl;
	cout<< " L = " << L << endl;
	cout<< " T = " << T << endl;
	cout<< " dt = " << dt <<endl;
	cout<< " rho = " << rho <<endl;
	cout<< " N = " << N <<endl;
	cout<< "\n" << stars <<endl;
	
	cout<< "Exact Solution:\n" <<endl;
	cout<< " J = rho(1-rho) = "<< rho*(1-rho) <<endl;
	cout<< " Density = rho = "<< rho <<endl;
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamD("");
    	fileNameStreamD << "R_Densprof" << "_rho="<< rho <<".txt";
    	std::string fileNameD = fileNameStreamD.str();
    	ofstream q1(fileNameD.c_str());
	
	#pragma omp parallel
	#pragma omp for lastprivate(lattice, J)

	initialise();
	update();

	for (j = 1; j <= L; j++) // Time averaging of the density.
	{
		SiteDens = densprof[j]/Tdif;
		q1<<j<<" "<< SiteDens <<endl;
		AvgDens += SiteDens;
	}
	
	q1.close();
        cputime2 = cputime ();
        cputime0 = cputime2 - cputime1;
	cout<< "Simulation Results:\n" <<endl;
	cout<<" J = "<< J/Tdif <<endl;
	cout<<" Density = "<< AvgDens/L <<endl;
	cout<< "\n";
	cout<< " Elapsed cpu time for main computation:\n";
	cout<< "  " << cputime2 << " seconds" <<endl;
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamL("");
    	fileNameStreamL << "R_Log" << "_rho="<< rho <<".txt";
    	std::string fileNameL = fileNameStreamL.str();
    	ofstream q2(fileNameL.c_str());
	
	q2<< stars <<endl;
	q2<< "Parameters:\n" <<endl;
	q2<< " L = " << L << endl;
	q2<< " T = " << T << endl;
	q2<< " dt = " << dt <<endl;
	q2<< " rho = " << rho <<endl;
	q2<< " N = " << N <<endl;
	q2<< stars <<endl;
	
	q2<< "Exact Solution:\n" <<endl;
	q2<< " J = rho(1-rho) = "<< rho*(1-rho) <<endl;
	q2<< " Density = rho = "<< rho <<endl;
	q2<< "\n" << stars <<endl;
	
	q2<< "Simulation Results:\n" <<endl;
	q2<<" J = "<< J/Tdif <<endl;
	q2<<" Density = "<<AvgDens/L<<endl;
	q2<< "\n";
	q2<< " Elapsed cpu time for main computation:\n";
	q2<< "  " << cputime2 << " seconds" <<endl;
	q2<< "\n" << stars <<endl;
	q2.close();
}
 
