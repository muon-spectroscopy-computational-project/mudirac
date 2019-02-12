/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * integrate.hpp
 * 
 * Routines for integrating differential equations with the shooting method
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <math.h>
#include <vector>
#include "utils.hpp"
#include "constants.hpp"

using namespace std;

void shoot_QP(vector<double> &Q, vector<double> &P, vector<double> AA, vector<double> AB, vector<double> BA, vector<double> BB,
              double h = 1, int stop_i = -1, char dir = 'f');