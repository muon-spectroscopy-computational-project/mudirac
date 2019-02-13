/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * boundary.cpp
 * 
 * Boundary conditions
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "boundary.hpp"

void boundaryDiracCoulomb(vector<double> &Q, vector<double> &P, vector<double> r, double E, int k, double m)
{
    int l = k < 0 ? -k - 1 : k; // Angular momentum number
    int N = Q.size();

    if (P.size() != N || r.size() != N || N < 4) {
        throw "Invalid array size passed to boundaryDiracCoulomb";
    }
    
}