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

CoulombSpherePotential::CoulombSpherePotential(double Z, double R)
{
    this->Z = Z;
    this->R = R;
    R3 = pow(R, 3);
    this->VR = R > 0 ? -1.5 * Z / R : 0;
}

double CoulombSpherePotential::V(double r)
{
    if (r < 0)
    {
        throw invalid_argument("Negative radius not allowed for CoulombPotential");
    }
    else if (r < R)
    {
        return Z * pow(r, 2) / (2 * R3) + VR;
    }
    else
    {
        return -Z / r;
    }
}

CoulombUehlingSpherePotential::CoulombUehlingSpherePotential(double Z, double R, int usteps) : CoulombSpherePotential(Z, R)
{
    this->usteps = usteps;
}

double CoulombUehlingSpherePotential::ukernel_r_greater(double u, double r, double R)
{
    double ans = exp(-2 * r * Physical::c / u);
    ans *= (exp(2 * R * Physical::c / u) * (R * u * Physical::alpha / 2 - pow(u * Physical::alpha, 2) / 4) +
            exp(-2 * R * Physical::c / u) * (R * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4));

    return ans;
}

double CoulombUehlingSpherePotential::ukernel_r_smaller(double u, double r, double R)
{
    double ans = (exp(-2 * r * Physical::c / u) - exp(2 * r * Physical::c / u));
    ans *= (exp(-2 * R * Physical::c / u) * (R * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4) -
            pow(u * Physical::alpha, 2) / 4);

    return ans;
}
