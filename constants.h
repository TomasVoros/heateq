#ifndef constants
#define constants

#include <vector>

// stefan-boltzmann
const double sigma = 0;//5.67e-8;
const double T0 = 280;                     // initial temperature

// meteoroid material parameters
const double Tb = 1000;                    // temperature T at beginning height
const double lam0 = lamt(T0);              // at 280K
const double c0 = 450;                     // at 280K
const double delta0 = 7800;                // at 280K
const double alpha0 = alphat(T0);          // at 280K

// meteoroid kinematics and atmosphere
const double Lam = 1;
const double vinf = 30000;
const double rho0 = 1.3;
const double H = 7000;
const double zr = 0;

// fitted constants for thermal diffusivity
//           brd[materials][borders]
//           fit[materials][segments][constants]
const int    abrd[1][1] = {{1068}};
const double afit[1][2][3] = {{{-2.156e-5, 4.342e-2, 682.2}, {6.61e-6, -8.844e-5, -1041}}};
const int    lbrd[1][2] = {{479, 1127}};
const double lfit[1][3][2] = {{{113.8, -0.1116}, {84.82, -0.05}, {16.76, 9.924e-3}}};

#endif
