#include "planet.h"

planet::planet() // Constructors of planets
{
	mass = 1.;
	position[0] = 1.;
	position[1] = 0.;
	velocity[0] = 0.;
	velocity[1] = 0.;
	potential = 0.;
	kinetic = 0.;
}

planet::planet(double M, double x, double y, double vx, double vy)
{
	mass = M;
	position[0] = x;
	position[1] = y;
	velocity[0] = vx;
	velocity[1] = vy;
	potential = 0.;
	kinetic = 0.;
}


double planet::distance(planet otherPlanet) // Returns distance between two planets
{
	double x1, y1, x2, y2, xx, yy;

	x1 = this->position[0];
	y1 = this->position[1];

	x2 = otherPlanet.position[0];
	y2 = otherPlanet.position[1];

	xx = x1 - x2;
	yy = y1 - y2;

	return sqrt(xx*xx + yy * yy);
}

double planet::KineticEnergy() //Returns kinetic energy of a given planet
{
	double velocity2 = (this->velocity[0] * this->velocity[0]) + (this->velocity[1] * this->velocity[1]);
	return 0.5*this->mass*velocity2;
}

double planet::PotentialEnergy(planet &otherPlanet, double Gconst) //Returns potential energy of a given planet
{
	return -Gconst * this->mass*otherPlanet.mass / this->distance(otherPlanet);
}