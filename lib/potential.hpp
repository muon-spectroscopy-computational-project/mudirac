/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * potential.hpp
 *
 * Classes for various types of potentials - header file
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
 */

#include <cmath>
#include <stdexcept>
#include <vector>

#include "constants.hpp"
#include "econfigs.hpp"
#include "integrate.hpp"
#include "utils.hpp"

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
  virtual double V(double r) { return 0; };
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
  virtual double V(double r) override;

protected:
  double R, R3, VR, Z;
};

/**
 * @brief Coulomb potential with 2-terms Fermi charge distribution
 * @note  Coulomb potential assuming a 2-terms Fermi charge distribution, using
 *        parameter c derived from R, and assuming a constant skin thickness t
 *        [L. Visscher, ATOMIC DATA AND NUCLEAR DATA TABLES 67, 207–224 (1997)]
 *
 * @retval None
 */
class CoulombFermi2Potential : public CoulombSpherePotential
{
public:
  CoulombFermi2Potential(double Z = 1.0, double R = -1, double A = 1.0,
                         double thickness = Physical::fermi2_T, int csteps = 5000);
  virtual double V(double r) override;

  double getc() { return c; }

protected:
  double c, innerV, dx, T;
  vector<vector<double>> grid;
  vector<double> Vgrid;
};

/**
 * @brief  Uehling correction term to Coulomb potential
 * @note   Uehling correction term to Coulomb potential (polarisation of the
 * vacuum, including only electron-positron pairs) assuming a uniform spherical
 * charge of radius R
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
  double ukernel_r_greater(int i, double r, bool rR = false);

  static double ukernel_r_smaller(double u, double r, double R);
  double ukernel_r_smaller(int i, double r, bool rR = false);

  static double ukernel_r_verysmall(double u, double R);
  static double ukernel_point(double u, double r);

protected:
  double exp_cutoff_high =
      INFINITY; // Cutoff point x beyond which we approximate exp(-x) = 0
  double exp_cutoff_low =
      0.0; // Cutoff point x under which we approximate exp(-x) = 1
  double Z, R, rho, K, V0, du, uint0;
  int usteps;
  vector<double> uarg, uker, u24c2, uker_great, uker_small;
};

/**
 * @brief  Coulomb potential from a charge distribution on a logarithmic grid
 * @note   Coulomb potential from a charge distribution on a fixed logarithmic
 * grid, interpolated for values that are not part of it.
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
  pair<int, int> getGridLimits() { return pair<int, int>(i0, i1); };

protected:
  int i0, i1;
  double rc, dx, Q, rho0, V0;
  vector<double> rho, Vpot;
  vector<vector<double>> grid;

  void initPotential(vector<double> rho);
};

/**
 * @brief  Coulomb potential from an electronic distribution on a logarithmic
 * grid
 * @note   Coulomb potential from a charge distribution on a fixed logarithmic
 * grid generated from an electronic configuration. Changes only the way it's
 * initialised, the potential is then calculated as in BkgGridPotential
 *
 * @retval None
 */
class EConfPotential : public BkgGridPotential
{
public:
  EConfPotential(ElectronicConfiguration econf = ElectronicConfiguration(),
                 double rc = 1, double dx = 1e-2, double rho_eps = 1e-5,
                 double max_r0 = -1, double min_r1 = -1);

protected:
  ElectronicConfiguration ec;
};

#endif