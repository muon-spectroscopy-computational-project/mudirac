/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * wavefunction.hpp
 *
 * Functions useful for operations on wavefunctions - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <math.h>
#include <vector>
#include "utils.hpp"
#include "integrate.hpp"

using namespace std;

double braKetLog(vector<double> psiL, vector<double> psiR, vector<double> r, double dx);
double braKetLog(vector<vector<double>> psiL, vector<vector<double>> psiR, vector<double> r, double dx);

double braOpKetLog(vector<double> psiL, vector<double> op, vector<double> psiR, vector<double> r, double dx);
double braOpKetLog(vector<double> psiL, vector<double> (*op)(vector<double> psi), vector<double> psiR, vector<double> r, double dx);

double braOpKetLog(vector<vector<double>> psiL, vector<double> op, vector<vector<double>> psiR, vector<double> r, double dx);
double braOpKetLog(vector<vector<double>> psiL, vector<vector<double>> (*op)(vector<vector<double>> psi),
                   vector<vector<double>> psiR, vector<double> r, double dx);