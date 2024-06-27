/**
 * MuDirac - A muonic atom Dirac equation solver
 * by Simone Sturniolo (2019-2020)
 *
 * potential.hpp
 *
 * Classes for various types of potentials
 *
 * @author Simone Sturniolo
 * @version 1.0 20/03/2020
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
CoulombSpherePotential::CoulombSpherePotential(double Z, double R) {

  LOG(INFO) << "Initialising Coulomb sphere potential, Z = " << Z
            << ", R = " << R << "\n";

  this->Z = Z;
  this->R = R;
  R3 = pow(R, 3);
  this->VR = R > 0 ? -1.5 * Z / R : 0;
}

double CoulombSpherePotential::V(double r) {

  if (r < 0) {
    throw invalid_argument("Negative radius not allowed for CoulombPotential");
  } else if (r < R) {
    return Z * pow(r, 2) / (2 * R3) + VR;
  } else {
    return -Z / r;
  }
}

/**
 * @brief  Initialise Coulomb Fermi 2-term potential
 * @note   Initialise Coulomb potential with
 * nuclear charge distributed according to a Fermi 2-term function
 *
 * @param  Z:         Nuclear charge
 * @param  R:         Nuclear radius
 * @param  A:         Atomic mass
 * @param  csteps:    Steps used to integrate numerically the potential
 * @param  thickness: Thickness parameter of the 'skin' of the nucleus
 * @retval
 */
CoulombFermi2Potential::CoulombFermi2Potential(double Z, double R, double A,
    double thickness, int csteps): CoulombSpherePotential(Z, R) {

  vector<double> rho;

  T = thickness;

  if (R < 0) {
    // Skip any initialisation, but the potential won't work properly
    return;
  }

  // First, define C for this radius
  if (A >= 5.0) {
    c = sqrt(R * R -
             7.0 / 3.0 * pow(M_PI * T / (4 * log(3.0)), 2));
  } else {
    c = 2.2291e-5 * pow(A, 1.0 / 3.0) - 0.90676e-5;
  }

  LOG(INFO) << "Initialising Coulomb Fermi-2 potential, c = " << c << "\n";

  // Then find the grid
  grid = logGrid(1e-8, 1e-2, csteps);
  // And define the density over it
  for (int i = 0; i < csteps; ++i) {
    rho.push_back(
      4.0 * M_PI * pow(grid[1][i], 2.0) /
      (1.0 + exp(4 * log(3.0) * (grid[1][i] - c) / T)));
  }
  // Integrate to find the total charge
  dx = grid[0][1] - grid[0][0];
  double Q = trapzInt(dx, vectorOperation(grid[1], rho, '*'));
  Q += 1.0 / 3.0 * grid[1][0] * rho[0]; // Inner sphere volume
  // Adjust to make sure that total charge is indeed Z
  rho = vectorOperation(rho, Z / Q, '*');
  innerV = rho[0] / (6 * pow(grid[1][0], 2.0));

  // Integrate potential
  Vgrid = vector<double>(csteps);
  shootPotentialLog(Vgrid, rho, dx);

  // Adjust potential so that it connects well at grid edge
  VR = -Z / grid[1].back() - Vgrid.back();
}

double CoulombFermi2Potential::V(double r) {

  if (r < 0) {
    throw invalid_argument("Negative radius not allowed for CoulombPotential");
  } else if (r > grid[1].back()) {
    return -Z / r;
  } else if (r < grid[1][0]) {
    return innerV * pow(r, 2) + VR;
  } else {
    double x = log(r / grid[1][0]);
    double xi = x / dx;
    int i0 = floor(xi);
    int i1 = ceil(xi);

    if (i0 == i1) {
      return Vgrid[i0] + VR;
    } else {
      // Interpolate
      return lerp(Vgrid[i0], Vgrid[i1], xi - i0) + VR;
    }
  }
  return 0.0;
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
UehlingSpherePotential::UehlingSpherePotential(double Z, double R, int usteps) {
  this->Z = Z;
  this->R = R;
  this->usteps = usteps;
  du = 1.0 / (usteps - 1.0);
  uarg = vector<double>(usteps, 0);
  uker = vector<double>(usteps, 0);
  u24c2 = vector<double>(usteps, 0);
  uker_great = vector<double>(usteps, 0);
  uker_small = vector<double>(usteps, 0);

  // Initialise values to be reused
  for (int i = 1; i < usteps; ++i) {
    double u = i * du;
    uker[i] = sqrt(1 - u * u) * (1 + 0.5 * u * u);
    u24c2[i] = pow(u * Physical::alpha, 2.0) / 4.0;
    uker_great[i] = ukernel_r_greater(u, 0, R);
    uker_small[i] = exp(-2 * R * Physical::c / u) *
                    (R * u * Physical::alpha / 2 + u24c2[i]) -
                    u24c2[i];
  }

  if (R > 0) {
    rho = Z * 0.75 / (M_PI * pow(R, 3));
    // Compute the 'uint0' term
    uarg[0] = 0;
    double eps = 0.5 * 1e-7 * du;
    for (int i = 1; i < usteps; ++i) {
      double u = i * du;
      uarg[i] = ukernel_r_verysmall(u, R) * uker[i];
    }
    uint0 = trapzInt(du, uarg);
  } else {
    rho = Z /
          (M_PI *
           Physical::alpha); // Works so that the front constant stays the same
  }

  K = -2 * pow(Physical::alpha, 2) / 3 * rho;
}

/**
 * @brief  Compute the Uehling integral kernel for r > R
 *
 * @param  u
 * @param  r
 * @param  R
 * @retval Kernel value
 */
double UehlingSpherePotential::ukernel_r_greater(double u, double r, double R) {
  double ans = exp(-2 * r * Physical::c / u);
  ans *= (exp(2 * R * Physical::c / u) *
          (R * u * Physical::alpha / 2 - pow(u * Physical::alpha, 2) / 4) +
          exp(-2 * R * Physical::c / u) *
          (R * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4));

  return ans;
}

/**
 * @brief  Compute the Uehling integral kernel for r > R
 * @note   Compute the Uehling integral kernel for r > R.
 * This version of the function makes use of stored, pre-computed
 * values for the u-grid and simplifications for the r = R case
 * for maximum efficiency.
 *
 * @param  i    Grid index for u (u = du*i)
 * @param  r
 * @param  rR   If true, consider r = R
 * @retval Kernel value
 */
double UehlingSpherePotential::ukernel_r_greater(int i, double r, bool rR) {
  double ans;
  if (!rR) {
    ans = exp(-2 * r * Physical::c / (du * i));
    ans *= uker_great[i];
  } else {
    double u = du * i;
    ans = (r * u * Physical::alpha / 2 - u24c2[i]) +
          exp(-4 * r * Physical::c / u) *
          (r * u * Physical::alpha / 2 + u24c2[i]);
  }
  return ans;
}

/**
 * @brief  Compute the Uehling integral kernel for r < R
 *
 * @param  u
 * @param  r
 * @param  R
 * @retval Kernel value
 */
double UehlingSpherePotential::ukernel_r_smaller(double u, double r, double R) {
  double ans = (exp(-2 * r * Physical::c / u) - exp(2 * r * Physical::c / u));
  ans *= (exp(-2 * R * Physical::c / u) *
          (R * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4) -
          pow(u * Physical::alpha, 2) / 4);

  return ans;
}

/**
 * @brief  Compute the Uehling integral kernel for r < R
 * @note   Compute the Uehling integral kernel for r < R.
 * This version of the function makes use of stored, pre-computed
 * values for the u-grid and simplifications for the r = R case
 * for maximum efficiency.
 *
 * @param  i    Grid index for u (u = du*i)
 * @param  r
 * @param  rR   If true, consider r = R
 * @retval Kernel value
 */
double UehlingSpherePotential::ukernel_r_smaller(int i, double r, bool rR) {
  double ans;
  double u = du * i;

  if (!rR) {
    ans = (exp(-2 * r * Physical::c / u) - exp(2 * r * Physical::c / u));
    ans *= uker_small[i];
  } else {
    double ecru = exp(-2 * r * Physical::c / u);
    ans = (r * u * Physical::alpha / 2 + u24c2[i]) * (pow(ecru, 2) - 1) +
          u24c2[i] * (1 / ecru - ecru);
  }
  return ans;
}

double UehlingSpherePotential::ukernel_r_verysmall(double u, double R) {
  return 4 * Physical::c / u *
         (-exp(-2 * R * Physical::c / u) *
          (0.5 * R * u * Physical::alpha + pow(u / (2 * Physical::c), 2)) +
          pow(u / (2 * Physical::c), 2));
}

double UehlingSpherePotential::ukernel_point(double u, double r) {
  return 1 / u * exp(-2 * r * Physical::c / u);
}

double UehlingSpherePotential::V(double r) {
  // Avoid all this mess if r is big enough
  if (r > exp_cutoff_high * 0.5 * Physical::alpha) {
    return 0.0;
  } else if (r <= exp_cutoff_low * 0.5 * du * Physical::alpha) {
    return K * uint0;
  }
  // Fill in the u integration kernel
  uarg[0] = 0;
  for (int i = 1; i < usteps; ++i) {
    double u = i * du;
    if (R <= 0) {
      uarg[i] = ukernel_point(u, r);
    } else if (r > R) {
      uarg[i] = ukernel_r_greater(i, r);
    } else {
      uarg[i] = ukernel_r_greater(i, r, true) + ukernel_r_smaller(i, r) -
                ukernel_r_smaller(i, r, true);
    }

    uarg[i] *= uker[i];
  }
  double ans = trapzInt(du, uarg);

  return K / r * ans;
}

BkgGridPotential::BkgGridPotential(vector<double> rho, double rc, double dx,
                                   int i0, int i1) {
  this->rc = rc;
  this->dx = dx;
  this->i0 = i0;
  this->i1 = i1;

  grid = logGrid(rc, dx, i0, i1);

  initPotential(rho);
}

BkgGridPotential::BkgGridPotential() {
  rc = 1;
  dx = 1e-3;
  i0 = 0;
  i1 = 0;
  Vpot = vector<double>(1, 0.0);
  V0 = 0;
  Q = 0;
  rho0 = 0;
  grid = logGrid(rc, dx, i0, i1);
}

void BkgGridPotential::initPotential(vector<double> rho) {
  this->rho = rho;
  rho0 = rho[0];
  Vpot = vector<double>(i1 - i0 + 1);
  shootPotentialLog(Vpot, rho, dx);
  // The charge is an integral plus an assumed constant charge density for r <
  // r0
  Q = trapzInt(dx, vectorOperation(rho, grid[1], '*')) +
      rho0 * grid[1][0] / 3.0;
  V0 = -Q / grid[1][i1 - i0] - Vpot[i1 - i0];
}

double BkgGridPotential::V(double r) {
  // Find the index
  double xi = log(r / rc) / dx;

  if (xi < i0) {
    return 1.0 / 6.0 * rho0 * pow(r / grid[1][0], 2.0) + V0;
  } else if (xi > i1) {
    return -Q / r;
  } else {
    // Interpolate
    int il = floor(xi);
    int ir = ceil(xi);

    if (il == ir) {
      return Vpot[il - i0] + V0;
    }
    double drl = grid[1][il - i0] * (exp((xi - il) * dx) - 1);
    double f = drl / (grid[1][ir - i0] - grid[1][il - i0]);
    return lerp(Vpot[il - i0], Vpot[ir - i0], f) + V0;
  }
}

double BkgGridPotential::Vgrid(int i) {
  if (i >= i0 && i <= i1) {
    return Vpot[i - i0] + V0;
  } else {
    double r = rc * exp(i * dx);
    if (i < i0) {
      return 1.0 / 6.0 * rho0 * pow(r / grid[1][0], 2.0) + V0;
    } else {
      return -Q / r;
    }
  }
}

/**
 * @brief  Initialise electronic configuration potential
 * @note   Initialise electronic configuration potential from an
 * ElectronicConfiguration object, grid parameters, and a tolerance. The limits
 * of the grid will then be established based on where the charge goes below the
 * given tolerance.
 *
 * @param  econf:       The electronic configuration to use
 * @param  rc:          Central point of the grid
 * @param  dx:          Logarithmic step
 * @param  rho_eps:     Tolerance value of the density at which to stop
 * acquiring it
 * @param  max_r0:      Maximum value for the inner radius of the grid. Will be
 * ignored if negative
 * @param  min_r1:      Minimum value for the outer radius of the grid. Will be
 * ignored if negative
 * @retval
 */
EConfPotential::EConfPotential(ElectronicConfiguration econf, double rc,
                               double dx, double rho_eps, double max_r0,
                               double min_r1) {
  this->ec = econf;
  this->rc = rc;
  this->dx = dx;

  // Now to find the boundaries...
  double r = rc;
  vector<double> rho = vector<double>(1, ec.hydrogenicChargeDensity(rc));
  i0 = 0;
  i1 = 0;

  max_r0 = max_r0 < 0 ? 2 * rc : max_r0;
  min_r1 = min_r1 < 0 ? rc / 2 : min_r1;

  while (abs(rho.front()) > rho_eps || r > max_r0) {
    i0--;
    r = rc * exp(i0 * dx);
    rho.insert(rho.begin(), ec.hydrogenicChargeDensity(r));
  }
  while (abs(rho.back()) > rho_eps || r < min_r1) {
    i1++;
    r = rc * exp(i1 * dx);
    rho.push_back(ec.hydrogenicChargeDensity(r));
  }

  LOG(INFO) << "Electronic configuration potential grid boundaries found:\n";
  LOG(INFO) << " i0 = " << i0 << " = " << rc * exp(i0 * dx) << "\n";
  LOG(INFO) << " i1 = " << i1 << " = " << rc * exp(i1 * dx) << "\n";

  grid = logGrid(rc, dx, i0, i1);

  // Now compute the charge density
  initPotential(rho);
}