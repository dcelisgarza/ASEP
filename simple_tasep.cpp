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
//================================================================================//
// Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP) //
//================================================================================//

// For details on the exact solution see:
// B. Derrida, M.R. Evans, V. Hakim, V. Pasquier, Exact solution of a 1d asymmetric exclusion model using a matrix formulation J. Phys. A26, 1493-1517 (1993) 

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
long double Teq = T*0.2;
long double Tdif = T - Teq;
long double dt = 0.01; 
long double t;

long int site;
long int nextsite;

long double J;
long double SiteDens;
long double AvgDens;

long double alpha;
long double beta;
long double densprof[L+1], lattice[L+1], p[L+1];
long double enter, eject, hop;

long int j;

//////////////////////////////////////// Functions ///////////////////////////////////////
void initialise() // Initialize lattice.
{
	J = 0.0; // Initial current
	AvgDens = 0.0;
	for (j = 1; j <= L; j++)
	{
		lattice[j]=0; densprof[j]=0;  p[j]=1;
	}
	enter = dt * alpha * p[1]; // The probability of entering.
	eject = dt * beta * p[L]; // The probability of exiting.
	hop = dt * p[L/2]; // The probability of hopping.
}
void move()
{	
	for (j = 1;j <= L; j++) // Goes through as many iterations as the array's length.
	{
		site = rand()%L + 1; // Picks out a random site in the array.
		nextsite = site +1; // Defines the next site. Saves the program from unecessary calculations.
		if (site == 1 && lattice[site] == 0 && ran3(&dum) <= enter) // Enter first site.
		{
			lattice[site] = 1; // If the criteria are met, the site is filled.
		}
		if (site == L && lattice[site] == 1 && ran3(&dum) <= eject) // Exit last site.
		{
			lattice[site] = 0; // If the citeria are met, the site is emptied.
		}		
		if (lattice[site] == 1 && lattice[nextsite] == 0 && site >=1 && site <= L-1 && ran3(&dum) <= hop) // Bulk of the array.
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
	for (t = 0; t <= T; t += dt)
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
	const string program ="Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP).";
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
	
        // display paramteres 	
	cout<< "\n" << stars <<endl;
	cout<< "Parameters:\n" <<endl;
	cout<< " L = " << L << endl;
	cout<< " T = " << T << endl;
	cout<< " dt = " << dt <<endl;
	cout<< " alpha = " << alpha << endl;
	cout<< " beta = " << beta << endl;
	cout<< "\n" << stars <<endl;
	
	cout<< "Exact Solution:\n" <<endl;
	cout<< "B. Derrida, M.R. Evans, V. Hakim, V. Pasquier, Exact solution of a 1d asymmetric exclusion model using a matrix formulation J. Phys. A26, 1493-1517 (1993)" <<endl;
	if (alpha > 0.5 && beta > 0.5)
	{
		cout<< " Maximal Current Phase" <<endl;
		cout<< "\n";
		cout<< " J = 0.25" <<endl;
		cout<< " Bulk Density = 0.5" <<endl;
	}
	if (beta < alpha && beta < 0.5)
	{
		cout<< " High Density Phase" <<endl;
		cout<< "\n";
		cout<< " J = b(1-b) = " << beta*(1-beta) <<endl;
		cout<< " Bulk Density = " << 1.0-beta <<endl;
	}
	if (alpha < beta && alpha < 0.5)
	{
		cout<< " Low Density Phase"<<endl;
		cout<< "\n";
		cout<< " J = a(1-a) = " << alpha*(1-alpha) <<endl;
		cout<< " Bulk Density = " << alpha <<endl;
	}
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamD("");
    	fileNameStreamD << "Densprof" <<"_a=" << alpha << "_b=" << beta << ".txt";
    	std::string fileNameD = fileNameStreamD.str();
    	ofstream q1(fileNameD.c_str());
	
	#pragma omp parallel
	#pragma omp for lastprivate(lattice, J)

	initialise();
	update();

	for (j = 1; j <= L; j++) // Time averaging of the density. It sweeps through the array.
	{
		SiteDens = densprof[j]/Tdif;
		q1<< j <<" "<< SiteDens <<endl;
		AvgDens += SiteDens;
	}
	
	q1.close();
        cputime2 = cputime ();
        cputime0 = cputime2 - cputime1;
	cout<< "Simulation Results:\n" <<endl;
	cout<< " J = "<< J/Tdif <<endl;
	cout<< " Bulk Density = "<< AvgDens/L <<endl;
	cout<< "\n";
	cout<< " Elapsed cpu time for main computation:\n";
	cout<< "  " << cputime2 << " seconds" <<endl;
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamL("");
    	fileNameStreamL << "Log" << "_a=" << alpha << "_b=" << beta <<".txt";
    	std::string fileNameL = fileNameStreamL.str();
    	ofstream q2(fileNameL.c_str());
	
	q2<< stars <<endl;
	q2<< "Parameters:\n" <<endl;
	q2<< " L = " << L << endl;
	q2<< " T = " << T << endl;
	q2<< " dt = " << dt <<endl;
	q2<< " alpha = " << alpha << endl;
	q2<< " beta = " << beta << endl;
	q2<< "\n" << stars << endl;
	
	q2<< "Exact Solution:\n" <<endl;
	q2<< "B. Derrida, M.R. Evans, V. Hakim, V. Pasquier, Exact solution of a 1d asymmetric exclusion model using a matrix formulation J. Phys. A26, 1493-1517 (1993)" <<endl;
	if (alpha > 0.5 && beta > 0.5)
	{
		q2<< " Maximal Current Phase" <<endl;
		q2<< "\n";
		q2<< " J = 0.25" <<endl;
		q2<< " Bulk Density = 0.5" <<endl;
	}
	if (beta < alpha && beta < 0.5)
	{
		q2<< " High Density Phase" <<endl;
		q2<< "\n";
		q2<< " J = b(1-b) = " << beta*(1-beta) <<endl;
		q2<< " Bulk Density = " << 1.0-beta <<endl;
	}
	if (alpha < beta && alpha < 0.5)
	{
		q2<< " Low Density Phase"<<endl;
		q2<< "\n";
		q2<< " J = a(1-a) = " << alpha*(1-alpha) <<endl;
		q2<< " Bulk Density = " << alpha <<endl;
	}
	q2<< "\n" << stars <<endl;
	
	q2<< "Simulation Results:\n" <<endl;
        q2<< " J = "<< J/Tdif <<endl;
	q2<< " Bulk Density = "<< AvgDens/L <<endl;
	q2<< "\n";
	q2<< " Elapsed cpu time for main computation:\n";
	q2<< "  " << cputime2 << " seconds" <<endl;
	q2<< "\n" << stars <<endl;
	q2.close();
}
 
