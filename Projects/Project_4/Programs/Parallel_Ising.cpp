#include "mpi.h"
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


// Inline function for Periodic boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
	return (i + limit + add) % (limit);
}
// Function to initialize tha lattice and its energy and magnetization
void InitializeLattice(int, mat &, double&, double&);
// The Monte Carlo method with the Metropolis algorithm 
void MetropolisSampling(int, int, double, vec &);
// Function that prints to file the results of the calculations  
void WriteResultstoFile(int, int, double, vec);

// Main program. Execute with "outputfile #spins #MCcycles InitialTemp FinalTemp TempStep" as arguments.

int main(int argc, char* argv[])
{
	string filename;
	int NSpins, MonteCarloCycles;
	double InitialTemp, FinalTemp, TempStep;
	int NProcesses, RankProcess;
	cout << "Estamos trabajando en ello." << endl;
	//  MPI initializations
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NProcesses);
	MPI_Comm_rank(MPI_COMM_WORLD, &RankProcess);
	if (RankProcess == 0 && argc <= 5) {
		cout << "Bad Usage: " << argv[0] <<
			" Please include output file, Number of spins, MC cycles, initial and final temperature and tempurate step as arguments" << endl;
		exit(1);
	}
	if ((RankProcess == 0) && (argc > 1)) {
		filename = argv[1];
		NSpins = atoi(argv[2]);
		MonteCarloCycles = atoi(argv[3]);
		InitialTemp = atof(argv[4]);
		FinalTemp = atof(argv[5]);
		TempStep = atof(argv[6]);
	}

	// Declare new file name and add lattice size to file name, only master node opens file
	if (RankProcess == 0) {
		string fileout = filename;
		string argument = to_string(NSpins);
		fileout.append(argument);
		ofile.open(fileout);
	}
	// Broadcast to all nodes common variables
	MPI_Bcast(&MonteCarloCycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&NSpins, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&TempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Start Monte Carlo sampling by looping over the temperature steps
	double  TimeStart, TimeEnd, TotalTime;
	TimeStart = MPI_Wtime(); //Start computing time
	for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature += TempStep) {
		vec LocalExpectationValues = zeros<mat>(5);
		// Start Monte Carlo computation and get local expectation values
		MetropolisSampling(NSpins, MonteCarloCycles, Temperature, LocalExpectationValues);
		// Create an array to allocate our results
		vec TotalExpectationValues = zeros<mat>(5);
		for (int i = 0; i < 5; i++) {
			MPI_Reduce(&LocalExpectationValues[i], &TotalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (RankProcess == 0) WriteResultstoFile(NSpins, MonteCarloCycles*NProcesses, Temperature, TotalExpectationValues);
	}
	if (RankProcess == 0)  ofile.close();  // close output file
	TimeEnd = MPI_Wtime(); //End of time computing
	TotalTime = TimeEnd - TimeStart;
	if (RankProcess == 0) {
		cout << "Time = " << TotalTime << " on number of processors: " << NProcesses << endl;
	}
	// End MPI
	MPI_Finalize();
	return 0;
}


// The Monte Carlo method with the Metropolis algorithm
void MetropolisSampling(int NSpins, int MonteCarloCycles, double Temperature, vec &ExpectationValues)
{
	// Initialize the seed and call the Mersienne algorithm
	std::random_device rd;
	std::mt19937_64 gen(rd());
	// Set up the uniform distribution for x \in [[0, 1]
	std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);
	// Initialize the lattice spin values
	mat SpinMatrix = zeros<mat>(NSpins, NSpins);
	// Initialize energy and magnetization 
	double Energy = 0.;     double MagneticMoment = 0.;
	// Initialize array for expectation values
	InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
	// Set up an array that includes the possible energy changes
	vec EnergyDifference = zeros<mat>(17);
	for (int de = -8; de <= 8; de += 4) EnergyDifference(de + 8) = exp(-de / Temperature);
	// Start Monte Carlo cycles
	int AllSpins = NSpins * NSpins;
	for (int cycles = 1; cycles <= MonteCarloCycles; cycles++) {
		// The sweep over the lattice, looping over all spin sites
		for (int Spins = 0; Spins < AllSpins; Spins++) {
			int ix = (int)(RandomNumberGenerator(gen)*NSpins);
			int iy = (int)(RandomNumberGenerator(gen)*NSpins);
			int deltaE = 2 * SpinMatrix(ix, iy)*
				(SpinMatrix(ix, PeriodicBoundary(iy, NSpins, -1)) +
					SpinMatrix(PeriodicBoundary(ix, NSpins, -1), iy) +
					SpinMatrix(ix, PeriodicBoundary(iy, NSpins, 1)) +
					SpinMatrix(PeriodicBoundary(ix, NSpins, 1), iy));
			if (RandomNumberGenerator(gen) <= EnergyDifference(deltaE + 8)) { //Metropolis algorithm acceptance condition
				SpinMatrix(ix, iy) *= -1.0;  // flip one spin and accept new spin configuration
				MagneticMoment += 2.0*SpinMatrix(ix, iy);
				Energy += (double)deltaE;
			}
		}
		// Update expectation values  for local node after a sweep through the lattice
		ExpectationValues(0) += Energy;    ExpectationValues(1) += Energy * Energy;
		ExpectationValues(2) += MagneticMoment;
		ExpectationValues(3) += MagneticMoment * MagneticMoment;
		ExpectationValues(4) += fabs(MagneticMoment);
	}
} // End of Metropolis sampling over spins

// Function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix, double& Energy, double& MagneticMoment)
{
	// Set up spin matrix and initial magnetization with all spins pointing up or down
	for (int x = 0; x < NSpins; x++) {
		for (int y = 0; y < NSpins; y++) {
			SpinMatrix(x, y) = 1.0; // Spin orientation for the ground state
			MagneticMoment += (double)SpinMatrix(x, y);
		}
	}
	// Set up initial energy
	for (int x = 0; x < NSpins; x++) {
		for (int y = 0; y < NSpins; y++) {
			Energy -= (double)SpinMatrix(x, y)*
				(SpinMatrix(PeriodicBoundary(x, NSpins, -1), y) +
					SpinMatrix(x, PeriodicBoundary(y, NSpins, -1)));
		}
	}
}// End of the initializing function



void WriteResultstoFile(int NSpins, int MonteCarloCycles, double temperature, vec ExpectationValues)
{
	double norm = 1.0 / ((double)(MonteCarloCycles));  // divided by  number of cycles 
	double E_ExpectationValues = ExpectationValues(0)*norm;
	double E2_ExpectationValues = ExpectationValues(1)*norm;
	double M_ExpectationValues = ExpectationValues(2)*norm;
	double M2_ExpectationValues = ExpectationValues(3)*norm;
	double Mabs_ExpectationValues = ExpectationValues(4)*norm;
	// all expectation values are per spin, divide by 1/NSpins/NSpins
	double AllSpins = 1.0 / ((double)NSpins*NSpins);
	double HeatCapacity = (E2_ExpectationValues - E_ExpectationValues * E_ExpectationValues)*AllSpins / temperature / temperature;
	double MagneticSusceptibility = (M2_ExpectationValues - Mabs_ExpectationValues * Mabs_ExpectationValues)*AllSpins / temperature;
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setw(15) << setprecision(8) << temperature;
	ofile << setw(15) << setprecision(8) << E_ExpectationValues * AllSpins;
	ofile << setw(15) << setprecision(8) << HeatCapacity;
	ofile << setw(15) << setprecision(8) << M_ExpectationValues * AllSpins;
	ofile << setw(15) << setprecision(8) << MagneticSusceptibility;
	ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues * AllSpins << endl;
} // End of the output function
