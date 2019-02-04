/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * utils.hpp
 * 
 * Utility functions of various sorts - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <vector>

using namespace std;

double effectiveMass(double m1, double m2);

int factorial(int n);

double genLaguerrePoly(double x, int k, int alpha);
vector<double> genLaguerrePoly(vector<double> x, int k, int alpha);