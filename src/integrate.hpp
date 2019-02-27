/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * integrate.hpp
 * 
 * Routines for integrating differential equations with the shooting method - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <math.h>
#include <vector>
#include "utils.hpp"
#include "constants.hpp"

using namespace std;

double trapzInt(vector<double> x, vector<double> y);
void shootQ(vector<double> &Q, vector<double> A, vector<double> B, double h = 1, int stop_i = -1, char dir = 'f');
void shootQP(vector<double> &Q, vector<double> &P, vector<double> AA, vector<double> AB, vector<double> BA, vector<double> BB,
             double h = 1, int stop_i = -1, char dir = 'f');
void shootNumerov(vector<double> &Q, vector<double> A, vector<double> B, double h = 1, int stop_i = -1, char dir = 'f');
void shootPotentialLog(vector<double> &V, vector<double> rho, double h = 1);

struct TurningPoint
{
    int i;
    double Qi, Qe, Pi, Pe;
};

TurningPoint shootDiracLog(vector<double> &Q, vector<double> &P, vector<double> r, vector<double> V,
                           double E, int k = -1, double m = 1, double dx = 1);
void shootDiracErrorDELog(vector<double> &zeta, vector<double> y, vector<double> r, vector<double> V,
                          int turn_i, double E, int k = -1, double m = 1, double dx = 1, char dir = 'f');