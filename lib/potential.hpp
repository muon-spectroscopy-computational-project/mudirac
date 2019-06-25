/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019)
 * 
 * potential.hpp
 * 
 * Classes for various types of potentials - header file
 * 
 * @author Simone Sturniolo
 * @version 0.1 04/02/2019
 */

#include <cmath>
#include <vector>
#include <stdexcept>
#include "integrate.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

#ifndef MUDIRAC_POTENTIAL
#define MUDIRAC_POTENTIAL

/**
  * @brief  Generic Potential class other potentials inherit from 
  * @note   A generic Potential class for other specific potentials
  * to inherit from.
  * 
  * @retval None
 */
class Potential
{
public:
  /**
    * @brief  Evaluate potential at r
    * @note   Evaluate potential at r
    * 
    * @param  r:  Distance from the center at which to evaluate potential
    * @retval     Potential
   */
  virtual double V(double r){};
};

/**
  * @brief  Coulomb potential with uniform spherical charge distribution
  * @note   Coulomb potential assuming a uniform spherical charge of radius R
  * 
  * @retval None
  */
class CoulombSpherePotential : Potential
{
public:
  CoulombSpherePotential(double Z = 1.0, double R = -1);
  double V(double r) override;

protected:
  double R, R3, VR, Z;
};

/**
  * @brief  Uehling correction term to Coulomb potential
  * @note   Uehling correction term to Coulomb potential (polarisation of the vacuum, 
  * including only electron-positron pairs) assuming a uniform spherical charge of radius R 
  *  
  * @retval None
 */
class UehlingSpherePotential : Potential
{
public:
  UehlingSpherePotential(double Z = 1.0, double R = -1, int usteps = 100);
  double V(double r) override;

  static double ukernel_r_greater(double u, double r, double R);
  static double ukernel_r_smaller(double u, double r, double R);
  static double ukernel_point(double u, double r);

protected:
  double Z, R, rho, K;
  int usteps;
  vector<double> uarg;
};

/**
  * @brief  Coulomb potential from a charge distribution on a logarithmic grid
  * @note   Coulomb potential from a charge distribution on a fixed logarithmic grid,
  * interpolated for values that are not part of it.
  * 
  * @retval None
 */
class BkgGridPotential : Potential
{
public:
  BkgGridPotential(vector<double> rho, double rc, double dx, int i0, int i1);
  double V(double r) override;
  double Vgrid(int i);

  double getQ() { return Q; };

protected:
  int i0, i1;
  double rc, dx, Q, rho0, V0;
  vector<double> Vpot;
  vector<vector<double>> grid;
};

#endif