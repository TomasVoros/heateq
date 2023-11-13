#ifndef functions
#define functions

#include <vector>

using namespace std;

// atmosphere density in height h
double rho(double);

// atmosphere density in height h, logarithmically interpolated
double rhox(double, const vector <double> &);

// heat flux in height h
double j(double, const vector <double> &);

// thermal diffusivity at temperature T
double alphat(double);

// heat conductivity at temperature T
double lamt(double);

// integral of heat conductivity
double t2theta(double);

// inverse of the previous function
double theta2t(double);

// reads atmospheric density from density.dat from 50km to 250km (including) with 1km step
void readdensity(vector <double> &);

#endif
