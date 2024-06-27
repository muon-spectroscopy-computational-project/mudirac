/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * utils.hpp
 *
 * Utility functions of various sorts - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
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
string printIupacState(int n, int l, bool s);
void parseIupacRange(string irange, vector<int> &nrange, vector<int> &lrange, vector<bool> &srange);

vector<string> splitString(string s, string sep = " ", bool merge = false, int maxn = -1);
string stripString(string s, string strip = " \t\n");
string upperString(string s);

// Functions useful for debugging
void writeTabulated2ColFile(vector<double> col1, vector<double> col2, string fname);

// Templated function declarations
/**
 * @brief  Search if a vector contains an element
 *
 * @tparam T        Type of the vectors' elements
 * @param  v        Vector
 * @param  e        Element
 * @retval          True if e is in v
 */
template <typename T>
bool vectorContains(vector<T> v, T e) {
  return find(v.begin(), v.end(), e) != v.end();
}

/**
 * @brief Perform an element-wise operation on two vectors
 *
 * @tparam T        Type of the vectors' elements
 * @param v1        First vector
 * @param v2        Second vector
 * @param op        Operator (can be +, -, *, / or ^)
 * @retval vector<T>
 */
template <typename T>
vector<T> vectorOperation(vector<T> v1, vector<T> v2, char op) {
  if (v1.size() != v2.size()) {
    throw runtime_error("Vectors don't match in size");
  }

  vector<T> ans(v1.size());

  T(*opfunc)
  (T, T);

  switch (op) {
    case '+':
      opfunc = [](T a, T b) {
        return a + b;
      };
      break;
    case '-':
      opfunc = [](T a, T b) {
        return a - b;
      };
      break;
    case '*':
      opfunc = [](T a, T b) {
        return a * b;
      };
      break;
    case '/':
      opfunc = [](T a, T b) {
        return a / b;
      };
      break;
    case '^':
      opfunc = pow;
      break;
    default:
      throw invalid_argument("Invalid operator code for vectorOperation");
      break;
  }

  for (int i = 0; i < v1.size(); ++i) {
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
 * @retval vector<T>
 */
template <typename T>
vector<T> vectorOperation(vector<T> v, T x, char op) {

  vector<T> ans(v.size());

  T(*opfunc)
  (T, T);

  switch (op) {
    case '+':
      opfunc = [](T a, T b) {
        return a + b;
      };
      break;
    case '-':
      opfunc = [](T a, T b) {
        return a - b;
      };
      break;
    case '*':
      opfunc = [](T a, T b) {
        return a * b;
      };
      break;
    case '/':
      opfunc = [](T a, T b) {
        return a / b;
      };
      break;
    case '^':
      opfunc = pow;
      break;
    default:
      throw invalid_argument("Invalid operator code for vectorOperation");
      break;
  }

  for (int i = 0; i < v.size(); ++i) {
    ans[i] = opfunc(v[i], x);
  }

  return ans;
}

#endif