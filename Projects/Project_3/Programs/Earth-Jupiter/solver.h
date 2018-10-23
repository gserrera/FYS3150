#ifndef SOLVER_H
#define SOLVER_H
#include "planet.h"
#include <vector>
#include <fstream>
using std::vector;

class solver
{
public:
	friend class planet;

	// properties
	double total_mass, G;
	int total_planets;
	vector<planet> all_planets;
	double totalKinetic;
	double totalPotential;

	// constants

	// initializers
	solver();

	// functions
	void add(planet newplanet);
	void print_position(std::ofstream &output, int dimension, double time, int number);
	void print_energy(std::ofstream &output, double time);
	void VelocityVerlet(int dimension, int integration_points, double final_time, int print_number);
	double **setup_acc(int height, int width);
	void delete_acc(double **matrix);
	void GravitationalForce(planet &current, planet &other, double &Fx, double &Fy);
	void KineticEnergySystem();
	void PotentialEnergySystem();
};

#endif // SOLVER_H