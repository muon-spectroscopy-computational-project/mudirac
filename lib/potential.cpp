/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * potential.hpp
 * 
 * Classes for various types of potentials
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include "potential.hpp"

CoulombPotential::CoulombPotential(double Z, double R)
{
    this->Z = Z;
    this->R = R;
    R3 = pow(R, 3);
    this->VR = R > 0 ? -1.5 * Z / R : 0;
}

double CoulombPotential::V(double r)
{
    if (r < 0)
    {
        throw invalid_argument("Negative radius not allowed for CoulombPotential");
    }
    else if (r < R)
    {
        return Z * pow(r, 2) / (2 * R3) + VR;
    }
    else {
        return -Z/r;
    }
}

