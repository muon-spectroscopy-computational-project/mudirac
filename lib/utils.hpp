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
#include <stdexcept>
#include <functional>

#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_UTILS
#define MUDIRAC_UTILS

double effectiveMass(double m1, double m2);

int factorial(int n);
double sinc(double x);
double lerp(double a, double b, double t);

template <typename T>
vector<T> vectorOperation(vector<T> v1, vector<T> v2, char op);
template <typename T>
vector<T> vectorOperation(vector<T> v, T x, char op);

vector<double> linGrid(double x0, double x1, int n = 100);
vector<vector<double>> logGrid(double x0, double x1, int n = 100);
vector<vector<double>> logGrid(double xc, double dx, int i0, int i1);

double genLaguerrePoly(double x, int n, double alpha);
vector<double> genLaguerrePoly(vector<double> x, int n, int alpha);

int countNodes(vector<double> v, double tol = 1e-8);

void qnumSchro2Dirac(int l, bool s, int &k);
void qnumDirac2Schro(int k, int &l, bool &s);
void qnumNodes2Principal(int nodes, int l, int &n);
void qnumPrincipal2Nodes(int n, int l, int &nodes);

double cgCoeff(int k, double m, bool s);

void parseIupacState(string istate, int &n, int &l, bool &s);

vector<string> splitString(string s, string sep = " ", bool merge = false);
string stripString(string s, string strip = " \t\n");

// Templated function declarations
/**
 * @brief Perform an element-wise operation on two vectors
 * 
 * @tparam T        Type of the vectors' elements
 * @param v1        First vector
 * @param v2        Second vector
 * @param op        Operator (can be +, -, *, / or ^)
 * @return vector<T> 
 */
template <typename T>
vector<T> vectorOperation(vector<T> v1, vector<T> v2, char op)
{
    if (v1.size() != v2.size())
    {
        throw runtime_error("Vectors don't match in size");
    }

    vector<T> ans(v1.size());

    T(*opfunc)
    (T, T);

    switch (op)
    {
    case '+':
        opfunc = [](T a, T b) { return a + b; };
        break;
    case '-':
        opfunc = [](T a, T b) { return a - b; };
        break;
    case '*':
        opfunc = [](T a, T b) { return a * b; };
        break;
    case '/':
        opfunc = [](T a, T b) { return a / b; };
        break;
    case '^':
        opfunc = pow;
        break;
    default:
        throw invalid_argument("Invalid operator code for vectorOperation");
        break;
    }

    for (int i = 0; i < v1.size(); ++i)
    {
        ans[i] = opfunc(v1[i], v2[i]);
    }

    return ans;
}

/**
 * @brief Perform an element-wise operation on a vector
 * and a scalar
 * 
 * @tparam T        Type of the vectors' elements
 * @param v         Vector
 * @param x         Scalar
 * @param op        Operator (can be +, -, *, / or ^)
 * @return vector<T> 
 */
template <typename T>
vector<T> vectorOperation(vector<T> v, T x, char op)
{

    vector<T> ans(v.size());

    T(*opfunc)
    (T, T);

    switch (op)
    {
    case '+':
        opfunc = [](T a, T b) { return a + b; };
        break;
    case '-':
        opfunc = [](T a, T b) { return a - b; };
        break;
    case '*':
        opfunc = [](T a, T b) { return a * b; };
        break;
    case '/':
        opfunc = [](T a, T b) { return a / b; };
        break;
    case '^':
        opfunc = pow;
        break;
    default:
        throw invalid_argument("Invalid operator code for vectorOperation");
        break;
    }

    for (int i = 0; i < v.size(); ++i)
    {
        ans[i] = opfunc(v[i], x);
    }

    return ans;
}

#endif