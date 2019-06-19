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

/**
  * @brief  Initialise Coulomb spherical potential
  * @note   Initialise Coulomb spherical potential with
  * uniform spherical nucleus
  * 
  * @param  Z:      Nuclear charge
  * @param  R:      Nuclear radius (if omitted or negative, point-like)
  * @retval 
  */
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

/**
  * @brief      Initialise Uehling potential with spherical charge
  * @note       Initialise Uehling potential correction to the Coulomb
  * potential, assuming a uniform spherical charge distribution at the
  * centre.
  * 
  * @param  Z:          Nuclear charge
  * @param  R:          Nuclear radius (if omitted or negative, point-like)
  * @param  usteps:     Steps used to integrate numerically over 'u' parameter
  * @retval 
  */
UehlingSpherePotential::UehlingSpherePotential(double Z, double R, int usteps)
{
    this->Z = Z;
    this->R = R;
    this->usteps = usteps;

    if (R > 0)
        rho = Z * 0.75 / (M_PI * pow(R, 3));
    else
        rho = Z / (M_PI * Physical::alpha); // Works so that the front constant stays the same

    K = -2 * pow(Physical::alpha, 2) / 3 * rho;

    uarg = vector<double>(usteps);
}

double UehlingSpherePotential::ukernel_r_greater(double u, double r, double R)
{
    double ans = exp(-2 * r * Physical::c / u);
    ans *= (exp(2 * R * Physical::c / u) * (R * u * Physical::alpha / 2 - pow(u * Physical::alpha, 2) / 4) +
            exp(-2 * R * Physical::c / u) * (R * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4));

    return ans;
}

double UehlingSpherePotential::ukernel_r_smaller(double u, double r, double R)
{
    double ans = (exp(-2 * r * Physical::c / u) - exp(2 * r * Physical::c / u));
    ans *= (exp(-2 * R * Physical::c / u) * (R * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4) -
            pow(u * Physical::alpha, 2) / 4);

    return ans;
}

double UehlingSpherePotential::ukernel_point(double u, double r)
{
    return 1 / u * exp(-2 * r * Physical::c / u);
}

double UehlingSpherePotential::V(double r)
{
    double du = 1.0 / (usteps - 1.0);
    // Fill in the u integration kernel
    uarg[0] = 0;
    for (int i = 1; i < usteps; ++i)
    {
        double u = i * du;
        if (R <= 0)
        {
            uarg[i] = ukernel_point(u, r);
        }
        else if (r > R)
        {
            uarg[i] = ukernel_r_greater(u, r, R);
        }
        else
        {
            uarg[i] = ukernel_r_greater(u, r, r) + ukernel_r_smaller(u, r, R) - ukernel_r_smaller(u, r, r);
        }
        
        uarg[i] *= sqrt(1 - u * u) * (1 + 0.5 * u * u);
    }
    double ans = trapzInt(du, uarg);

    return K / r * ans;
}
