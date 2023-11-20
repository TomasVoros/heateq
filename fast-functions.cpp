#include <cmath>
#include <vector>
#include <fstream>
#include "functions.h"
#include "constants.h"

using namespace std;

// atmosphere density in height h
double rho(double h)
{
    return rho0 * exp(-h/H);
}

// atmosphere density in height h, logarithmically interpolated
double rhox(double h, const vector <double> &density)
{
    if(h < 50000) return 1;                // if height is below 50km

    int a = floor(h / 1000) - 50;          // because density.dat has 201 values of density between 50km and 250km
    int b =  ceil(h / 1000) - 50;
    if(a == b) return density[a];          // if height in km is integer

    double fa = density[a];                // logarithmically interpolate
    double fb = density[b];
    double expa = ceil(h / 1000) - h / 1000;
    double expb = h / 1000 - floor(h / 1000);
    return pow(fa, expa) * pow(fb, expb);
}

// heat flux in height h
double j(double h, const vector <double> &density)
{
    //return 0.5 * Lam * pow(vinf, 3) * rho0 * exp(-h/H);
    return 0.5 * Lam * pow(vinf, 3) * rhox(h, density);
}

// thermal diffusivity at temperature T
double alphat(double T)
{
    if(T < 1068) return -2.156e-5 + 4.342e-2/(T + 682.2);
    else         return 6.61e-6 - 8.844e-5/(T - 1041);
}

// heat conductivity at temperature T
double lamt(double T)
{
    if(     T < 479)  return 113.8 - 0.1116*T;
    else if(T < 1127) return 84.82 - 0.05*T;
    else              return 16.76 + 9.924e-3*T;
}

// integral of heat conductivity
double t2theta(double T)
{
    if(T < 479)
    {
        return (113.8*(T - T0) - 0.5*0.1116*(T*T - T0*T0))/lam0;
    }
    else if(T < 1127)
    {
        double C1 = 113.8*(479 - T0) - 0.5*0.1116*(479*479 - T0*T0);
        return (84.82*(T - 479) - 0.5*0.05*(T*T - 479*479) + C1)/lam0;
    }
    else
    {
        double C1 = 113.8*(479 - T0) - 0.5*0.1116*(479*479 - T0*T0);
        double C2 = 84.82*(1127 - 479) - 0.5*0.05*(1127*1127 - 479*479);
        return (16.76*(T - 1127) + 0.5*9.924e-3*(T*T - 1127*1127) + C1 + C2)/lam0;
    }
}

// inverse of the previous function
double theta2t(double theta)
{
    if(theta < t2theta(479))
    {
        double C0 = 113.8*T0 - 0.5*0.1116*T0*T0;
        return 113.8/0.1116 - sqrt(113.8*113.8 - 2*0.1116*(theta*lam0 + C0))/0.1116;
    }
    else if(theta < t2theta(1127))
    {
        double C1 = 113.8*(479 - T0) - 0.5*0.1116*(479*479 - T0*T0);
        double C2 = 84.82*479 - 0.5*0.05*479*479;
        return 84.82/0.05 - sqrt(84.82*84.82 - 2*0.05*(theta*lam0 + C2 - C1))/0.05;
    }
    else
    {
        double C3 = 113.8*(479 - T0) - 0.5*0.1116*(479*479 - T0*T0);
        double C4 = 84.82*(1127 - 479) - 0.5*0.05*(1127*1127 - 479*479);
        double C5 = 16.76*1127 + 0.5*9.924e-3*1127*1127;
        return -16.76/(9.924e-3) + sqrt(16.76*16.76 + 2*9.924e-3*(theta*lam0 + C5 - C3 - C4))/(9.924e-3);
    }
}

// reads atmospheric density from density.dat from 50km to 250km (including) with 1km step
void readdensity(vector <double> &density)
{
    double dummy;
    ifstream fin("density.dat");
    while(fin >> dummy) density.push_back(dummy);

}
