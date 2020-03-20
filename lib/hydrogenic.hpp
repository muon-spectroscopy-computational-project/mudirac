/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 * 
 * hydrogenic.hpp
 * 
 * Computes known hydrogenic atom solutions for reference - header file
 * 
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <math.h>
#include <vector>
#include "utils.hpp"
#include "constants.hpp"

using namespace std;

#ifndef MUDIRAC_HYDROGENIC
#define MUDIRAC_HYDROGENIC

double hydrogenicSchroEnergy(double Z = 1.0, double mu = 1.0, int n = 1);
double hydrogenicSchroWavefunction(double r, double Z = 1, double mu = 1, int n = 1, int l = 0);
vector<double> hydrogenicSchroWavefunction(vector<double> r, double Z = 1, double mu = 1, int n = 1, int l = 0);

double hydrogenicDiracEnergy(double Z = 1.0, double mu = 1.0, int n = 1, int k = -1, bool bind = false);
vector<double> hydrogenicDiracWavefunction(double r, double Z = 1, double mu = 1, int n = 1, int k = -1);
vector<vector<double>> hydrogenicDiracWavefunction(vector<double> r, double Z = 1, double mu = 1, int n = 1, int k = -1);

#endif