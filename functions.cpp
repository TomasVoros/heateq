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
    int m = 0;
    if(T < abrd[m][0]) return afit[m][0][0] + afit[m][0][1]/(T + afit[m][0][2]);
    else               return afit[m][1][0] + afit[m][1][1]/(T + afit[m][1][2]);
}

// heat conductivity at temperature T
double lamt(double T)
{
    int m = 0;
    if(     T < lbrd[m][0]) return lfit[m][0][0] + lfit[m][0][1]*T;
    else if(T < lbrd[m][1]) return lfit[m][1][0] + lfit[m][1][1]*T;
    else                    return lfit[m][2][0] + lfit[m][2][1]*T;
}

// integral of heat conductivity
double t2theta(double T)
{
    int m = 0;
    double a0 = lfit[m][0][0],
           b0 = lfit[m][0][1],
           a1 = lfit[m][1][0],
           b1 = lfit[m][1][1],
           a2 = lfit[m][2][0],
           b2 = lfit[m][2][1];
    double brd1 = lbrd[m][0],
           brd2 = lbrd[m][1];

    if(T < brd1)
    {
        return (a0*(T - T0) + 0.5*b0*(T*T - T0*T0))/lam0;
    }
    else if(T < brd2)
    {
        double C1 = a0*(brd1 - T0) + 0.5*b0*(brd1*brd1 - T0*T0);
        return (a1*(T - brd1) + 0.5*b1*(T*T - brd1*brd1) + C1)/lam0;
    }
    else
    {
        double C1 = a0*(brd1 - T0) + 0.5*b0*(brd1*brd1 - T0*T0);
        double C2 = a1*(brd2 - brd1) + 0.5*b1*(brd2*brd2 - brd1*brd1);
        return (a2*(T - brd2) + 0.5*b2*(T*T - brd2*brd2) + C1 + C2)/lam0;
    }
}

// inverse of the previous function
double theta2t(double theta)
{
    int m = 0;
    double a0 = lfit[m][0][0],
           b0 = lfit[m][0][1],
           a1 = lfit[m][1][0],
           b1 = lfit[m][1][1],
           a2 = lfit[m][2][0],
           b2 = lfit[m][2][1];
    double brd1 = lbrd[m][0],
           brd2 = lbrd[m][1];

    if(theta < t2theta(brd1))
    {
        double C0 = a0*T0 + 0.5*b0*T0*T0;
        return -a0/b0 + sqrt(a0*a0 + 2*b0*(theta*lam0 + C0))/b0;
    }
    else if(theta < t2theta(brd2))
    {
        double C1 = a0*(brd1 - T0) + 0.5*b0*(brd1*brd1 - T0*T0);
        double C2 = a1*brd1 + 0.5*b1*brd1*brd1;
        return -a1/b1 + sqrt(a1*a1 + 2*b1*(theta*lam0 + C2 - C1))/b1;
    }
    else
    {
        double C3 = a0*(brd1 - T0) + 0.5*b0*(brd1*brd1 - T0*T0);
        double C4 = a1*(brd2 - brd1) + 0.5*b1*(brd2*brd2 - brd1*brd1);
        double C5 = a2*brd2 + 0.5*b2*brd2*brd2;
        return -a2/b2 + sqrt(a2*a2 + 2*b2*(theta*lam0 + C5 - C3 - C4))/b2;
    }
}

// reads atmospheric density from density.dat from 50km to 250km (including) with 1km step
void readdensity(vector <double> &density)
{
    double dummy;
    ifstream fin("density.dat");
    while(fin >> dummy) density.push_back(dummy);

}
