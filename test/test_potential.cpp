#include "../lib/atom.hpp"
#include "../lib/potential.hpp"
#include <iostream>
#include <math.h>
#include <vector>

#include "../vendor/catch/catch.hpp"

#define ERRTOL_HIGH 1e-5
#define ERRTOL_LOW 1e-3

TEST_CASE("Coulomb potential", "[CoulombSpherePotential]") {
  // Point-like
  CoulombSpherePotential cpot = CoulombSpherePotential(1);
  REQUIRE(cpot.V(0) == -INFINITY);
  // Finite size
  cpot = CoulombSpherePotential(1, 1.0);
  REQUIRE(cpot.V(1) == -1);
  REQUIRE(cpot.V(0) == Approx(-1.5));
  REQUIRE(cpot.V(2) == Approx(-0.5));
  REQUIRE_THROWS_AS(cpot.V(-1), invalid_argument);
}

TEST_CASE("Coulomb Fermi-2 potential", "[CoulombFermi2Potential]") {
  CoulombFermi2Potential cpot = CoulombFermi2Potential(1, 1e-5);

  // Start by checking c
  REQUIRE(cpot.getc() == Approx(1.32234e-5));

  // Far limit
  REQUIRE(cpot.V(2.0) == Approx(-0.5));         // r > R
  REQUIRE(cpot.V(9e-3) == Approx(-1.0 / 9e-3)); // c << r < R
  // Near limit
  REQUIRE(cpot.V(1e-9) == Approx(-40614.7308).epsilon(1e-3));

  // r ~ c
  REQUIRE(cpot.V(1e-5) == Approx(-38519.8597).epsilon(1e-3));
  REQUIRE(cpot.V(3.981e-5) == Approx(-23338.5869).epsilon(1e-3));

  // Inheritance test
  CoulombSpherePotential *point = new CoulombFermi2Potential(1, 1e-5);

  REQUIRE(point->V(1e-5) == Approx(-38519.8597).epsilon(1e-3));
  REQUIRE(point->V(3.981e-5) == Approx(-23338.5869).epsilon(1e-3));
}

TEST_CASE("Coulomb Fermi-2 potential, grid settings", "[CoulombFermi2Potential]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::warning,
                                      AixLog::Type::all);
  // Pb-208 with R = 7.102018661385413 fm, t = 2.3 fm
  double R = 7.102018661385413 * Physical::fm;
  double t = 2.3 * Physical::fm;

  // Reference values from data/fermi2_coulomb_ref.py (mpmath, 40 digits),
  // at r = x*c
  vector<pair<double, double>> refs = {
    {0.1, -939016.31956046096}, {0.5, -865170.92018901327},
    {1.0, -642145.94556580227}, {1.5, -435445.28385701454},
    {3.0, -217733.74136103255}
  };
  CoulombFermi2Potential vlin(82, R, 208, t, -1);
  CoulombFermi2Potential vherm(82, R, 208, t, -1, 5000, true);
  REQUIRE(vherm.getc() == vlin.getc());
  double maxerr_lin = 0;
  for (const pair<double, double> &ref : refs) {
    double r = ref.first * vlin.getc();
    // Gauss-Legendre moments with cubic Hermite read-back: near machine precision
    REQUIRE(fabs(vherm.V(r) / ref.second - 1) < 1e-11);
    maxerr_lin = max(maxerr_lin, fabs(vlin.V(r) / ref.second - 1));
  }
  // The default construction carries a grid error of order 1e-6
  REQUIRE(maxerr_lin > 1e-7);

  REQUIRE_THROWS_AS(CoulombFermi2Potential(82, R, 208, t, -1, 1), invalid_argument);

  // K1-L3 (2p3/2 - 1s1/2) of muonic Pb-208 through Atom::setFermi2Grid
  auto k1l3 = [&](int csteps, bool hermite) {
    DiracAtom da = DiracAtom(82, Physical::m_mu, 208, FERMI2, 7.102018661385413);
    da.setFermi2Grid(csteps, hermite);
    da.setFermi2(t);
    return (da.getState(2, 1, true).E - da.getState(1, 0, true).E) / Physical::eV;
  };
  double E_default = k1l3(5000, false);
  double E_lin_fine = k1l3(80000, false);
  double E_herm = k1l3(5000, true);
  double E_herm_fine = k1l3(80000, true);
  // The Hermite construction is converged at the default grid size...
  REQUIRE(fabs(E_herm - E_herm_fine) < 1e-3);
  // ...and agrees with the linear construction on a 16 times finer grid
  REQUIRE(fabs(E_herm - E_lin_fine) < 0.1);
  // while the default construction sits about 24 eV below (issue #86)
  REQUIRE(E_herm - E_default > 20.0);

  DiracAtom da = DiracAtom(82, Physical::m_mu, 208, FERMI2);
  REQUIRE_THROWS_AS(da.setFermi2Grid(1, false), invalid_argument);
}

TEST_CASE("Uehling correction to Coulomb potential",
          "[CoulombSpherePotential]") {
  UehlingSpherePotential cpot = UehlingSpherePotential(1, 1e-4, 100);

  // Test static functions
  REQUIRE(cpot.ukernel_r_greater(0.5, 1, 0) == Approx(0));
  REQUIRE(cpot.ukernel_r_smaller(0.5, 1, 0) == Approx(0));
  REQUIRE(cpot.ukernel_r_smaller(0.5, 0, 1) == Approx(0));
  REQUIRE(cpot.ukernel_point(0.5, 0) == Approx(2));

  // More complex cases
  REQUIRE(cpot.ukernel_r_greater(0.5, 0.5 / (2 * Physical::c),
                                 0.5 / (2 * Physical::c)) ==
          Approx(0.125 / pow(Physical::c * M_E, 2)));
  REQUIRE(cpot.ukernel_r_smaller(0.5, 0.5 / (2 * Physical::c),
                                 0.5 / (2 * Physical::c)) ==
          Approx(2 * sinh(1) * 0.125 / pow(Physical::c, 2) * (0.5 - exp(-1))));
  REQUIRE(cpot.ukernel_point(0.5, 0.5 / (2 * Physical::c)) ==
          Approx(2 * exp(-1)));

  // Now the actual potential
  // Point-like atom (U)
  cpot = UehlingSpherePotential(92, -1, 2000);
  REQUIRE(cpot.V(0.01) == Approx(-0.1795).epsilon(1e-3));
  REQUIRE(cpot.V(5e-5) == Approx(-10225.5698).epsilon(1e-3));
  // Finite size atom (U-238 with A^1/3 size model)
  cpot = UehlingSpherePotential(92, 1.2 * pow(238.0, 1.0 / 3.0) * Physical::fm,
                                500);
  REQUIRE(cpot.V(0.01) == Approx(-0.1796).epsilon(1e-3));
  REQUIRE(cpot.V(5e-5) == Approx(-4438.7471).epsilon(1e-3));

  // Now test that it still works with approximations
  double rlow = 0.5 * Physical::alpha * 1e-6;
  double Vlow = cpot.V(rlow);
  double rhigh = 0.5 * Physical::alpha * 10;
  double Vhigh = cpot.V(rhigh);
  cpot.set_exp_cutoffs(1e-3, 5.0);
  REQUIRE(cpot.V(rlow) == Approx(Vlow).epsilon(1e-3));
  REQUIRE(cpot.V(rhigh) ==
          Approx(Vhigh).margin(1e-3)); // This being a 0, we use absolute values
}

TEST_CASE("Background charge on grid potential", "[BkgGridPotential]") {
  double rc = 1, dx = 0.000001;
  int i0 = -100000, i1 = 10000;
  vector<vector<double>> grid = logGrid(rc, dx, i0, i1);
  vector<double> rho(grid[0].size());

  for (int i = 0; i < rho.size(); ++i) {
    rho[i] = 3 * pow(grid[1][i], 2) * (grid[1][i] <= 1);
  }

  BkgGridPotential bpot = BkgGridPotential(rho, rc, dx, i0, i1);

  REQUIRE(bpot.getQ() == Approx(1).epsilon(1e-3));

  // For constant density it should match the regular coulomb potential
  CoulombSpherePotential cpot = CoulombSpherePotential(1, 1);

  REQUIRE(bpot.Vgrid(i0) == Approx(cpot.V(grid[1][0])));
  REQUIRE(bpot.Vgrid(0) == Approx(cpot.V(grid[1][-i0])));
  REQUIRE(bpot.Vgrid(i1) == Approx(cpot.V(grid[1][i1 - i0])));

  // Check that the interpolation works well even near the grid points
  REQUIRE(bpot.Vgrid(i0) == Approx(bpot.V(grid[1][0])));
  REQUIRE(bpot.Vgrid(0) == Approx(bpot.V(grid[1][-i0])));
  REQUIRE(bpot.Vgrid(i1) == Approx(bpot.V(grid[1][i1 - i0])));

  // Check behaviour at long range
  REQUIRE(bpot.V(2) == Approx(-0.5));
}

TEST_CASE("Electronic configuration potential", "[EConfPotential]") {
  ElectronicConfiguration hydrogen("H");
  EConfPotential hpot(hydrogen, 1, 1e-2, 1e-3);

  vector<double> rgrid = hpot.getGrid()[1];

  REQUIRE(rgrid[0] < 2.3e-2);
  REQUIRE(rgrid[rgrid.size() - 1] > 5.5);

  // Now test what happens when you add in grid limits
  hpot = EConfPotential(hydrogen, 1, 1e-2, 1e-3, 1e-2, 6);
  rgrid = hpot.getGrid()[1];

  REQUIRE(rgrid[0] < 1e-2);
  REQUIRE(rgrid[rgrid.size() - 1] > 6);
}