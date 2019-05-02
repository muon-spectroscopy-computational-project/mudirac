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

#include <cmath>
#include <string>
#include <vector>

using namespace std;

#ifndef MUDIRAC_UTILS
#define MUDIRAC_UTILS

double effectiveMass(double m1, double m2);

int factorial(int n);
double sinc(double x);

vector<double> linGrid(double x0, double x1, int n = 100);
vector<vector<double>> logGrid(double x0, double x1, int n = 100);

double genLaguerrePoly(double x, int n, double alpha);
vector<double> genLaguerrePoly(vector<double> x, int n, int alpha);

int countNodes(vector<double> v, double tol = 1e-8);

vector<string> splitString(string s, string sep = " ", bool merge = false);
string stripString(string s, string strip = " \t\n");
#endif