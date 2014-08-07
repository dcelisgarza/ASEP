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
//==================================================================================//
// Monte Carlo Siumulation of partially asymmetric simple exclusion process (PASEP) //
//==================================================================================//

// For details on the exact solution see:
// R. A. Blythe and M. R. Evans. "Nonequilibrium steady state of matrix-product form: as solver's guide." J. Phys. A: Math. Theor. Vol. 40 (2007) R333-R441.

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
long double Tdif = T - Teq;
long double  dt = 0.01; 
long double  t;

int site;
int nextsite;
int prevsite;
int pprevsite;
int nnextsite;

long double J1;
long double J2;
long double SiteDens;
long double AvgDens;

long double alpha;
long double beta;
long double gama;
long double delta;
long double eta;
long double chi;
long double densprof[L+1], lattice[L+2], p[L+1], pn[L+1];
long double l_en, l_ex, r_en, r_ex, f_hop, b_hop;

int j;

////////////////////////Functions/////////////////////////
void initialise()
{
	J1 = 0.0;
	J2 = 0.0;
	AvgDens = 0.0;
	lattice[0] = 0;
	lattice[L+1] = 0;
	for(j = 1; j <= L; j++)
	{
		lattice[j]=0; densprof[j]=0;  p[j]=gama; pn[j]=delta;
	}
	l_en = dt*alpha*p[1];
	l_ex = dt*eta*pn[1];

	r_en = dt*chi*pn[L];
	r_ex = dt*beta*p[L];
	
	f_hop = dt*p[L/2];
	b_hop = dt*pn[L/2];	
}

void boundary_interactions()
{
	// Left
	if (site == 1)
	{
		// FORWARD
		if (lattice[site]==0 && ran3(&dum)<=l_en)
		{
			lattice[site]=1;
		} 
		//BACKWARD
		if (lattice[site]==1 && ran3(&dum)<=l_ex)
		{
			lattice[site]=0;
		}
	}
	// Right
	if (site == L)
	{
		// FORWARD
		if (lattice[site] == 1 && ran3(&dum) <= r_ex)
		{
			lattice[site] = 0;             // moves to state (1,0) 
		}
		// BACKWARD
		if (lattice[site] == 0 && ran3(&dum) <= r_en)
		{
			lattice[site] = 1;             // moves to state (1,1) 
		}
	}
}

void move()
{
	for(j = 1; j <= L; j++)
	{
		site = (rand() % L+1); 
		prevsite = site - 1;
		nextsite = site + 1;
		boundary_interactions();
		// Forward
		if (lattice[site] == 1 && lattice[nextsite] == 0 && site >= 1 && site <= L-1 && ran3(&dum) <= f_hop)
		{
			lattice[nextsite] = 1;      // moves to state (0,0,1,0) 
			lattice[site] = 0; 
			if (site == L/2 && t >= Teq)
			{	
				J1++;	
			}
		} // End of site availability check.
		//Backward
		if (lattice[site] == 1 && lattice[prevsite] == 0 && site >= 2 && site <= L && ran3(&dum) <= b_hop)
		{
			lattice[prevsite] = 1;      // moves to state (0,0,1,0) 
			lattice[site] = 0; 
			if (site == L/2 && t >= Teq)
			{
				J2++;	
			}
		} // End of site availability check.
	} // End of loop through array.
}

void update()
{
	dum=-time(NULL);
	ran3(&dum);
	for (t = 0.0; t <= T; t+=dt) // Time loop allows for better averaging, and allows the system to reach steady state.
	{
		move();
		if(t >= Teq)
		{
			for(j = 1;j <= L; j++)
			{
				if (lattice[j] == 1)
				{
					densprof[j] += dt;
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
	const string program ="Monte Carlo Siumulation of partially asymmetric simple exclusion process (PASEP).";
	const string spaces(program.size(), '*');
	const string stars = spaces;
	cout<<"\n"<<endl;
	cout<<stars<<endl;
	cout<< program <<endl;
	cout<<stars<<endl;
	
	cout<< "\nalpha = ";
	cin>> alpha;
	cin.ignore();
	
	cout<< "eta = ";
	cin>> eta;
	cin.ignore();

	cout<< "beta = ";
	cin>> beta;
	cin.ignore();
	
	cout<< "chi = ";
	cin>> chi;
	cin.ignore();

	cout<< "gamma = ";
	cin>> gama;
	cin.ignore();

	cout<< "delta = ";
	cin>> delta;
	cin.ignore();
	
	cout<< "\n" << stars <<endl;
	cout<< "Parameters:\n" <<endl;
	cout<< " L = " << L << endl;
	cout<< " T = " << T << endl;
	cout<< " dt = " << dt <<endl;
	cout<< " alpha = " << alpha << endl;
	cout<< " eta = " << eta << endl;
	cout<< " beta = " << beta << endl;
	cout<< " chi = " << chi << endl;
	cout<< " gamma = " << gama << endl;
	cout<< " delta = " << delta << endl;
	cout<< "\n" << stars <<endl;
	
	cout<< "Exact Solution:\n" <<endl;
	cout<< " R. A. Blythe and M. R. Evans." <<endl;
	cout<< " ''Nonequilibrium steady state of matrix-product form: as solver's guide.''" <<endl;
	cout<< " J. Phys. A: Math. Theor. Vol. 40 (2007) R333-R441.\n" <<endl;
	if (alpha > 0.5 && beta > (1-delta)/2.0)
	{
		cout<< " Maximal Current Phase" <<endl;
		cout<< "\n";
		cout<< " J = (1-d)/4 = " << (1-delta)/4.0 <<endl;
	}
	if (beta < alpha && beta < (1-delta)/2.0)
	{
		cout<< " High Density Phase" <<endl;
		cout<< "\n";
		cout<< " J = b(1-d-b)/(1-d) = " << beta*(1-delta-beta)/(1-delta) <<endl;
	}
	if (alpha < beta && alpha < (1-delta)/2.0)
	{
		cout<< " Low Density Phase"<<endl;
		cout<< "\n";
		cout<< " J = a(1-d-a)/(1-d) = " << alpha*(1-delta-alpha)/(1-delta) <<endl;
	}
	cout<< "\n" << stars <<endl;

	std::ostringstream fileNameStreamD("");
    	fileNameStreamD << "D" <<"_a" << alpha <<"et" << eta << "b" << beta << "c" << chi << "g" << gama << "d" << delta << ".txt";
    	std::string fileNameD = fileNameStreamD.str();
    	ofstream q1(fileNameD.c_str());
	
	#pragma omp parallel
	#pragma omp for lastprivate(lattice, J)
  
	initialise();
	update();

	for(j = 1; j <= L; j++)
	{
		SiteDens = densprof[j]/Tdif;
		q1<< j <<" "<< SiteDens <<endl;
		AvgDens += SiteDens;
	}
	q1.close();
        cputime2 = cputime ();
        cputime0 = cputime2 - cputime1;
	cout<< "Simulation Results:\n" <<endl;
	cout<< " J = "<< (J1-J2)/Tdif <<endl;
	cout<< " Bulk Density = "<< AvgDens/L <<endl;
	cout<< "\n";
	cout<< " Elapsed cpu time for main computation:\n";
	cout<< "  " << cputime2 << " seconds" <<endl;
	cout<< "\n" << stars <<endl;
	
	std::ostringstream fileNameStreamL("");
    	fileNameStreamL << "L" <<"_a" << alpha <<"et" << eta << "b" << beta << "c" << chi << "g" << gama << "d" << delta << ".txt";
    	std::string fileNameL = fileNameStreamL.str();
    	ofstream q2(fileNameL.c_str());
	
	q2<< stars <<endl;
	q2<< "Parameters:\n" <<endl;
	q2<< " L = " << L << endl;
	q2<< " T = " << T << endl;
	q2<< " dt = " << dt <<endl;
	q2<< " alpha = " << alpha << endl;
	q2<< " eta = " << eta << endl;
	q2<< " beta = " << beta << endl;
	q2<< " chi = " << chi << endl;
	q2<< " gamma = " << gama << endl;
	q2<< " delta = " << delta << endl;
	q2<< "\n" << stars << endl;
	
	q2<< "Exact Solution:\n" <<endl;
	q2<< " R. A. Blythe and M. R. Evans." <<endl;
	q2<< " ''Nonequilibrium steady state of matrix-product form: as solver's guide.''" <<endl;
	q2<< " J. Phys. A: Math. Theor. Vol. 40 (2007) R333-R441.\n" <<endl;
	if (alpha > 0.5 && beta > (1-delta)/2.0)
	{
		q2<< " Maximal Current Phase" <<endl;
		q2<< "\n";
		q2<< " J = (1-d)/4 = " << (1-delta)/4.0 <<endl;
	}
	if (beta < alpha && beta < (1-delta)/2.0)
	{
		q2<< " High Density Phase" <<endl;
		q2<< "\n";
		q2<< " J = b(1-d-b)/(1-d) = " << beta*(1-delta-beta)/(1-delta) <<endl;
	}
	if (alpha < beta && alpha < (1-delta)/2.0)
	{
		q2<< " Low Density Phase"<<endl;
		q2<< "\n";
		q2<< " J = a(1-d-a)/(1-d) = " << alpha*(1-delta-alpha)/(1-delta) <<endl;
	}
	q2<< "\n" << stars << endl;
	
	q2<< "Simulation Results:\n" <<endl;
	q2<< " J = "<< (J1-J2)/Tdif <<endl;
	q2<< " Bulk Density = "<< AvgDens/L <<endl;
	q2<< "\n";
	q2<< " Elapsed cpu time for main computation:\n";
	q2<< "  " << cputime2 << " seconds" <<endl;
	q2<< "\n" << stars << endl;
	q2.close();	
}
 
