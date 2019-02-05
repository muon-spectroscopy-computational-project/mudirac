/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * utils.cpp
 * 
 * Utility functions of various sorts
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "utils.hpp"

/**
 * @brief  Compute the effective mass of a two-body system
 * 
 * @param  m1:  First mass
 * @param  m2:  Second mass
 * @retval      Effective mass
 */
double effectiveMass(double m1, double m2) {
    return (m1*m2)/(m1+m2);
}

/**
 * @brief  Factorial of n
 * 
 * @param  n:   Argument
 * @retval      n!
 */
int factorial(int n) {
    int fn = 1;

    if (n < 0) {
        throw "Factorial of a negative number";
    }

    for (int i = 1; i <= n; ++i) {
        fn *= i;
    }

    return fn;
}

/**
 * @brief  Compute the value of a generalised Laguerre polynomial at x
 * 
 * @param  x:       Argument of the polynomial 
 * @param  n:       Parameter n of the polynomial
 * @param  alpha:   Parameter alpha of the polynomial
 * @retval          Value of the polynomial
 */
double genLaguerrePoly(double x, int n, double alpha) {
    // We do this recursively. Not the most efficient way, but since we only need small numbers, 
    // it should be sufficient.

    double L = 1.0;

    if (n < 0) {
        throw "Parameter k of generalised Laguerre polynomial must be >= 0.";
    } 

    switch (n) {
        case 0:
            return L;
            break;
        case 1:
            return (1.0 + alpha - x);
            break;
        default:
            return ((2*n-1+alpha-x)*genLaguerrePoly(x, n-1, alpha)-(n-1+alpha)*genLaguerrePoly(x, n-2, alpha))/n;
            break;
    }
}

/**
 * @brief  Compute the value of a generalised Laguerre polynomial at multiple x
 * 
 * @param  x:       Arguments of the polynomial 
 * @param  n:       Parameter n of the polynomial
 * @param  alpha:   Parameter alpha of the polynomial
 * @retval          Value of the polynomial
 */
vector<double> genLaguerrePoly(vector<double> x, int n, int alpha) {
    // Just iteratively go over all x
    vector<double> L(x.size());

    for (int i = 0; i < x.size(); ++i) {
        L[i] = genLaguerrePoly(x[i], n, alpha);
    }

    return L;
}
