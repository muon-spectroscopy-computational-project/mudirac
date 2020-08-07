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