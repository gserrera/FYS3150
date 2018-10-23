#include "solver.h"
#include "planet.h"
#include <iostream>
#include <cmath>
#include "time.h"

solver::solver() //Constructor
{
	total_planets = 0;
	total_mass = 0;
	G = 4 * M_PI*M_PI;
	totalKinetic = 0;
	totalPotential = 0;
}


void solver::add(planet newplanet) //Adds planet to the solver object
{
	total_planets += 1;
	total_mass += newplanet.mass;
	all_planets.push_back(newplanet);
}

void solver::print_position(std::ofstream &output, int dimension, double time, int number)
{   // Writes mass, position and velocity to a file "output"
	output << "Time" << "\t" << "Planet" << "\t" << "Mass" << "\t" << "X" << "\t" << "Y" << "\t" << "Vx" << "\t" << "Vy" << std::endl;
		for (int i = 0; i < number; i++) {
			planet &Current = all_planets[i];
			output << time << "\t" << i + 1 << "\t" << Current.mass;
			for (int j = 0; j < dimension; j++) output << "\t" << Current.position[j];
			for (int j = 0; j < dimension; j++) output << "\t" << Current.velocity[j];
			output << std::endl;
		
		}
}

void solver::print_energy(std::ofstream &output, double time)
{   // Writes energies to a file "output"
	output << "Time" << "\t" << "Planet" << "\t" << "Kenergy" << "\t" << "Penergy" << "\t" << "E" << "\t" << "FullsystE" << std::endl;
	this->KineticEnergySystem();
	this->PotentialEnergySystem();
	for (int nr = 0; nr < total_planets; nr++) {
		planet &Current = all_planets[nr];
		output << time << "\t" << nr << "\t";
		output << Current.kinetic << "\t" << Current.potential << "\t" << Current.kinetic + Current.potential << "\t" << totalKinetic + totalPotential << std::endl;
	}
}


void solver::Euler(int dimension, int integration_points, double final_time, int print_number)
{   //  Euler solver for two coupeled ODEs in a given number of dimensions.
	
	// Define time step
	double time_step = final_time / ((double)integration_points);
	double time = 0.0;

	// Create files for data storage
	char *filename = new char[1000];
	char *filenameE = new char[1000];
	sprintf(filename, "SolarSystemEu_%d_%.2f.txt", total_planets, time_step); //Appends the number of planets and the time step to the output file name
	sprintf(filenameE, "SolarSystemEnergyEu_%d_%.2f.txt", total_planets, time_step);
	std::ofstream output_file(filename);
	std::ofstream output_energy(filenameE);

	// Set up arrays
	double **acceleration = setup_acc(total_planets, 2);
	double **acceleration_new = setup_acc(total_planets, 2);

	// Initialize forces
	double Fx, Fy; // Forces in each dimension

	// Write initial values of position and energy to output file
	print_position(output_file, dimension, time, print_number);
	print_energy(output_energy, time);

	// Set up clock
	clock_t planet_Eu, finish_Eu;
	planet_Eu = clock();

	// Loop over time
	time += time_step;
	while (time < final_time) {

		// Loop over all planets
		for (int nr1 = 0; nr1 < total_planets; nr1++) {
			planet &current = all_planets[nr1]; // Current planet we are looking at

			Fx = Fy = 0.0; // Reset forces before each run

			// Calculate new position for current planet
			for (int j = 0; j < dimension; j++) {
				current.position[j] += current.velocity[j] * time_step;
			}

			// Loop over all other planets to update the forces
			for (int nr2 = 0; nr2 < total_planets; nr2++) {
				if (nr2 != nr1) {
					planet &other = all_planets[nr2];
					GravitationalForce(current, other, Fx, Fy);
				}
			}

			// Update the acceleration with current forces
			acceleration[nr1][0] = Fx / current.mass;
			acceleration[nr1][1] = Fy / current.mass;

			// Calculate new velocity for current planet
			for (int j = 0; j < dimension; j++) current.velocity[j] += time_step*acceleration[nr1][j];
		}

		// Write current values to file and increase time
		print_position(output_file, dimension, time, print_number);
		print_energy(output_energy, time);

		//Increases the simulation time
		time += time_step;
	}
	// Stop clock and print out time usage
	finish_Eu = clock();
	std::cout << "Total time = " << "\t" << ((float)(finish_Eu - planet_Eu) / CLOCKS_PER_SEC) << " seconds" << std::endl; // print elapsed time

	// Close files
	output_file.close();
	output_energy.close();

	// Clear memory
	delete_acc(acceleration);
	delete_acc(acceleration_new);
}

double ** solver::setup_acc(int height, int width)
{   // Function to set up the 2D acceleration array

	// Set up matrix
	double **matrix;
	matrix = new double*[height];

	// Allocate memory
	for (int i = 0; i < height; i++)
		matrix[i] = new double[width];

	// Set values to zero
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			matrix[i][j] = 0.0;
		}
	}
	return matrix;
}

void solver::delete_acc(double **matrix)
{   // Function to deallocate memory of the acceleration 2D array

	for (int i = 0; i < total_planets; i++)
		delete[] matrix[i];
	delete[] matrix;
}

void solver::GravitationalForce(planet &current, planet &other, double &Fx, double &Fy)
{   // Function that calculates the gravitational force between two objects, component by component.

	// Calculate relative distance in each component and the total distance between current planet and the other
	double relative_distance[2];

	for (int j = 0; j < 2; j++) relative_distance[j] = current.position[j] - other.position[j];
	double r = current.distance(other);


	// Calculate the forces in each direction
	Fx -= this->G*current.mass*other.mass*relative_distance[0] / (r*r*r);
	Fy -= this->G*current.mass*other.mass*relative_distance[1] / (r*r*r);
}

void solver::KineticEnergySystem() //Calculates the kinetic energy of both the current planet and the whole system
{
	totalKinetic = 0;
	for (int nr = 0; nr < total_planets; nr++) {
		planet &Current = all_planets[nr];
		Current.kinetic = Current.KineticEnergy();
		totalKinetic += Current.kinetic;
	}
}

void solver::PotentialEnergySystem() //Calculates the potential energy of both the current planet and the whole system
{
	totalPotential = 0;
	for (int nr = 0; nr < total_planets; nr++) {
		planet &Current = all_planets[nr];
		Current.potential = 0;
	}
	for (int nr1 = 0; nr1 < total_planets; nr1++) {
		planet &Current = all_planets[nr1];
		for (int nr2 = 0; nr2 < total_planets; nr2++) {
			if (nr2 != nr1) {
				planet &Other = all_planets[nr2];
				Current.potential += Current.PotentialEnergy(Other, G);
				totalPotential += Current.potential;
			}
		}
	}
}