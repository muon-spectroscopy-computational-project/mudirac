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
#include "econfigs.hpp"
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

  void set_exp_cutoffs(double low, double high)
  {
    exp_cutoff_low = low;
    exp_cutoff_high = high;
  };

  static double ukernel_r_greater(double u, double r, double R);
  static double ukernel_r_smaller(double u, double r, double R);
  static double ukernel_r_verysmall(double u, double R);
  static double ukernel_point(double u, double r);

protected:
  double exp_cutoff_high = INFINITY; // Cutoff point x beyond which we approximate exp(-x) = 0
  double exp_cutoff_low = 0.0;       // Cutoff point x under which we approximate exp(-x) = 1
  double Z, R, rho, K, V0, du, uint0;
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
  BkgGridPotential();
  BkgGridPotential(vector<double> rho, double rc, double dx, int i0, int i1);
  double V(double r) override;
  double Vgrid(int i);

  double getQ() { return Q; };
  double getV0() { return V0; };
  vector<double> getrho() { return rho; };
  vector<vector<double>> getGrid() { return grid; };

protected:
  int i0, i1;
  double rc, dx, Q, rho0, V0;
  vector<double> rho, Vpot;
  vector<vector<double>> grid;

  void initPotential(vector<double> rho);
};

/**
 * @brief  Coulomb potential from an electronic distribution on a logarithmic grid
 * @note   Coulomb potential from a charge distribution on a fixed logarithmic grid
 * generated from an electronic configuration. Changes only the way it's initialised,
 * the potential is then calculated as in BkgGridPotential
 * 
 * @retval None
 */
class EConfPotential : public BkgGridPotential
{
  public:
    EConfPotential(ElectronicConfiguration econf=ElectronicConfiguration(), double rc=1, double dx=1e-2, double rho_eps=1e-5, double max_r0=-1, double min_r1=-1);

  protected:
    ElectronicConfiguration ec;
};

#endif