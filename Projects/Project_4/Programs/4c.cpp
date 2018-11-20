#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

// inline function for Periodic boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
	return (i + limit + add) % (limit);
}
// Function to initialise the spin matrix and its energy and magnetization
void InitializeLattice(int, mat &, double&, double&, int);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &, int);
// Funtion that prints to file the results of the calculations  
void WriteResultstoFile(int, int, double, vec, int);

// Main program

int main(int argc, char* argv[])
{
	string filename;
	int NSpins, MCcycles, Random;
	double Temp;
	cout << "Introduce the name of the output file: ";
	cin >> filename;
	cout << "Introduce the number of spins n (nxn lattice): ";
	cin >> NSpins;
	cout << "Introduce the number of MC cycles: ";
	cin >> MCcycles;
	cout << "Introduce the temperature: ";
	cin >> Temp;
	cout << "Introduce the initial lattice (0 for oriented lattice, 1 for random): ";
	cin >> Random;

	// Declare new file name and add lattice size to file name
	string fileout = filename;
	string argument = to_string(NSpins);
	fileout.append(argument);
	ofile.open(fileout);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "       Cycle:           E:          VarE:       M:                accepted:" << endl;
	// Array that will allocate the expectation values
	vec ExpectationValues = zeros<mat>(3);
	// Start Monte Carlo computation and get expectation values
	MetropolisSampling(NSpins, MCcycles, Temp, ExpectationValues, Random);
	// 
	ofile.close();  // close output file
	cout << "Done";
	return 0;
}



// The Monte Carlo method which includes the Metropolis algorithm
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues, int Random)
{
	// Initialize the seed and call the Mersienne algo
	std::random_device rd;
	std::mt19937_64 gen(rd());
	// Set up the uniform distribution for x \in [[0, 1]
	std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);
	// Initialize the lattice spin values
	mat SpinMatrix = zeros<mat>(NSpins, NSpins);
	//    initialize energy and magnetization 
	double Energy = 0.;     double MagneticMoment = 0.;
	// initialize the lattice
	InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment, Random);
	// setup an array that includes the possible energy changes
	vec EnergyDifference = zeros<mat>(17);
	for (int de = -8; de <= 8; de += 4) EnergyDifference(de + 8) = exp(-de / Temperature);
	// Loop over Monte Carlo cycles
	for (int cycles = 1; cycles <= MCcycles; cycles++) {
		int accept = 0;
		// The sweep over the lattice, looping over all spin sites
		for (int x = 0; x < NSpins; x++) {
			for (int y = 0; y < NSpins; y++) {
				int ix = (int)(RandomNumberGenerator(gen)*(double)NSpins);
				int iy = (int)(RandomNumberGenerator(gen)*(double)NSpins);
				int deltaE = 2 * SpinMatrix(ix, iy)*
					(SpinMatrix(ix, PeriodicBoundary(iy, NSpins, -1)) +
						SpinMatrix(PeriodicBoundary(ix, NSpins, -1), iy) +
						SpinMatrix(ix, PeriodicBoundary(iy, NSpins, 1)) +
						SpinMatrix(PeriodicBoundary(ix, NSpins, 1), iy));
				if (RandomNumberGenerator(gen) <= EnergyDifference(deltaE + 8)) {
					SpinMatrix(ix, iy) *= -1.0;  // flip one spin and accept new spin config
					MagneticMoment += (double)2 * SpinMatrix(ix, iy);
					Energy += (double)deltaE;
					accept += 1;
				}
			}
		}
		// update expectation values  for local node

		ExpectationValues(0) += Energy;    
		ExpectationValues(1) += fabs(MagneticMoment);
		ExpectationValues(2) += Energy * Energy;
		
		WriteResultstoFile(NSpins, cycles, Temperature, ExpectationValues, accept);
	}
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix, double& Energy, double& MagneticMoment, int Random)
{
	if (Random == 0) { //Oriented initial lattice
		// setup spin matrix and initial magnetization
		for (int x = 0; x < NSpins; x++) {
			for (int y = 0; y < NSpins; y++) {
				SpinMatrix(x, y) = 1.0; // spin orientation for the ground state
				MagneticMoment += (double)SpinMatrix(x, y);
			}
		}
		// setup initial energy
		for (int x = 0; x < NSpins; x++) {
			for (int y = 0; y < NSpins; y++) {
				Energy -= (double)SpinMatrix(x, y)*
					(SpinMatrix(PeriodicBoundary(x, NSpins, -1), y) +
						SpinMatrix(x, PeriodicBoundary(y, NSpins, -1)));
			}
		}
	}
	else if (Random == 1) { //Random initial lattice

		// setup spin matrix and initial magnetization
		for (int x = 0; x < NSpins; x++) {
			for (int y = 0; y < NSpins; y++) {
				std::random_device rd;
				std::mt19937_64 gen(rd());
				std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);
				if (RandomNumberGenerator(gen) >= 0.5) {
					SpinMatrix(x, y) = 1.0; // spin orientation for the ground state
				}
				else {
					SpinMatrix(x, y) = -1.0; // spin orientation for the ground state
				}
				MagneticMoment += (double)SpinMatrix(x, y);
			}
		}
		// setup initial energy
		for (int x = 0; x < NSpins; x++) {
			for (int y = 0; y < NSpins; y++) {
				Energy -= (double)SpinMatrix(x, y)*
					(SpinMatrix(PeriodicBoundary(x, NSpins, -1), y) +
						SpinMatrix(x, PeriodicBoundary(y, NSpins, -1)));
			}
		}

	}
}// end of the initializing funtion



void WriteResultstoFile(int NSpins, int cycles, double temperature, vec ExpectationValues, int accept)
{
	double E_ExpectationValues = ExpectationValues(0)/cycles;
	double Mabs_ExpectationValues = ExpectationValues(1)/cycles;
	double E2_ExpectationValues = ExpectationValues(2) / cycles;
	double E_variance = (E2_ExpectationValues - E_ExpectationValues * E_ExpectationValues) / NSpins / NSpins;

	// all expectation values are per spin, divide by 1/NSpins/NSpins
	
	ofile << setw(15) << setprecision(8) << cycles;
	ofile << setw(15) << setprecision(8) << E_ExpectationValues / NSpins / NSpins;
	ofile << setw(15) << setprecision(8) << E_variance;
	ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues / NSpins / NSpins;
	ofile << setw(15) << setprecision(8) << accept << endl;
} // end output function