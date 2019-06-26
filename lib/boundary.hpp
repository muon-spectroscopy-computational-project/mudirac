/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * boundary.hpp
 * 
 * Boundary conditions - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <math.h>
#include <vector>
#include "constants.hpp"
#include "utils.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_BOUNDARY
#define MUDIRAC_BOUNDARY

void boundaryDiracCoulomb(vector<double> &Q, vector<double> &P, vector<double> r, double E, int k = -1, double m = 1, double Z = 1,
                          double R = -1);
void boundaryDiracErrorDECoulomb(vector<double> &zeta, double E, int k = -1, double m = 1);

#endif