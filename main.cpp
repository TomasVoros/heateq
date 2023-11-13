// coefficients are formally function of temperature with stefan-boltzmann, explicit method

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "heateq/functions.h"
#include "heateq/constants.h"

using namespace std;

int main()
{
    // simulation parameters
    double h = 2.5e5;                      // initial height
    const double l = 2000;                // space elements
    const double dx = 1/l;                 // space step
    const double dt = 1e-5;                // time step

    // arrays of temperature and diffusion coefficient in time t and t+dt
    vector <double> theta1(l, t2theta(T0));
    vector <double> theta2(l, t2theta(T0));
    vector <double> alpha(l, alpha0);

    // get atmospheric density
    vector <double> density;
    readdensity(density);

    int cnt = 0;
    while(theta2t(theta1[0]) < T0 + Tb)
    {
        h -= vinf*cos(zr)*dt;
        for(int i = 0; i < l; i++) alpha[i] = alphat(theta2t(theta1[i]));

        theta2[0] = theta1[0] + alpha[0]*dt/(dx*dx)*(theta1[1] - theta1[0]) + alpha[0]*dt/dx*(j(h, density) - sigma*(pow(theta2t(theta1[0]), 4) - pow(T0, 4)))/lam0;
        for(int i = 1; i < l-1; i++) theta2[i] = theta1[i] + alpha[i]*dt/(dx*dx)*(theta1[i+1] - 2*theta1[i] + theta1[i-1]);
        theta2[l-1] = theta1[l-1] + alpha[l-1]*dt/(dx*dx)*(theta1[l-2] - theta1[l-1]);

        for(int i = 0; i < l; i++) theta1[i] = theta2[i];
        cnt++;
        if(cnt % 667 == 0) cout << theta2t(theta1[0]) - T0 << '\t' << h/1000 << endl;
    }
    cout << theta2t(theta1[0]) - T0 << '\t' << h/1000 << endl;
    cout << "theory" << '\t' << -H*log(2*Tb/(Lam*rho0) * sqrt(lam0*delta0*c0*cos(zr)/(H*pow(vinf, 5)))) / 1000 << endl << endl;

    cout << "dt " << dt << endl;
    cout << "l " << l << endl;
    cout << "log LHS " << log10(2 * Tb * sqrt(lam0*c0*delta0) / Lam) << endl;
    cout << "log RHS " << log10(pow(vinf, 2.5) * rho(h) / sqrt(cos(zr)/H)) << endl;
    cout << "STAB " << 2*alpha0*dt/(dx*dx) << endl;
    cout << "last element " << theta2t(theta1[l-1]) << endl;

    //ofstream fout;
    //fout.open("data.dat");
    //for(int i = 0; i < l; i++) if(isinf(log10(T1[i] == 0))) fout << i << '\t' << log10(T1[i]) << endl;
    //for(int i = 280; i < 2000; i++) fout << i << '\t' << theta2t(t2theta(i)) << endl;
    //system("wgnuplot -persist -e \"plot 'data.dat' with lines\"");
    return 0;
}
