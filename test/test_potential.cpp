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

TEST_CASE("Uehling potential, cancellation-free finite-size kernel",
          "[UehlingSpherePotential]") {
  // Reference values from data/uehling_sphere_ref.py (mpmath, 40 digits)
  struct UehlingRef {
    int Z;
    double R_fm, x, V; // x = r/R
  };
  vector<UehlingRef> refs = {
    {82, 7.1, 0.01, -4420.8619722329687}, {82, 7.1, 0.5, -3992.0617712912353},
    {82, 7.1, 0.99, -2659.6613644582458}, {82, 7.1, 1.0, -2620.1596811696094},
    {82, 7.1, 1.5, -1447.0531678537094}, {82, 7.1, 4.0, -324.41078457627675},
    {14, 4.0, 0.01, -1581.234188791387}, {14, 4.0, 0.5, -1430.7353466425679},
    {14, 4.0, 0.99, -966.90696977011124}, {14, 4.0, 1.0, -953.29931105167961},
    {14, 4.0, 1.5, -542.97490468111887}, {14, 4.0, 4.0, -134.40890407494777}
  };

  for (const UehlingRef &ref : refs) {
    double R = ref.R_fm * Physical::fm;
    double r = ref.x * R;

    // Gauss-Legendre in theta reaches the reference with few nodes
    // (r close to R needs more nodes, so it gets a looser bound)
    UehlingSpherePotential ugauss(ref.Z, R, 400, false, true);
    double tol_gauss = (ref.x == 0.99) ? 5e-8 : 1e-10;
    REQUIRE(fabs(ugauss.V(r) / ref.V - 1) < tol_gauss);

    // Trapezoidal rule at many steps: limited only by the step size
    // (the original form of the kernel loses all precision inside the
    // nucleus here)
    UehlingSpherePotential ustable(ref.Z, R, 16000, true);
    REQUIRE(fabs(ustable.V(r) / ref.V - 1) < 1.5e-6);
  }

  // The trapezoidal error shrinks monotonically as the steps grow, also inside
  // the nucleus, where the original form of the kernel breaks down
  vector<int> steps = {1000, 1200, 1400, 2000, 4000, 8000, 16000};
  for (const UehlingRef &ref : refs) {
    if (ref.Z != 82 || ref.x < 0.5 || ref.x > 1.0)
      continue;
    double R = ref.R_fm * Physical::fm;
    double r = ref.x * R;
    double err_prev = INFINITY;
    for (int U : steps) {
      UehlingSpherePotential ustable(ref.Z, R, U, true);
      double err = ustable.V(r) / ref.V - 1;
      REQUIRE(err < 0);
      REQUIRE(fabs(err) < err_prev);
      err_prev = fabs(err);
    }
  }

  // The central value (r = 0, and below a low cutoff) uses the selected rule.
  // With 2 nodes the trapezoidal grid holds only u = 0 and u = 1, where the
  // weighted central kernel vanishes, while the Gauss rule does not: the value
  // must be the Gauss sum of the r -> 0 limit of the kernel divided by r.
  {
    double Rc = 7.1 * Physical::fm;
    double K = -2 * pow(Physical::alpha, 2) / 3 * (82 * 0.75 / (M_PI * pow(Rc, 3)));
    double V2 = 0;
    for (double xg : {-1 / sqrt(3.0), 1 / sqrt(3.0)}) {
      double th = (xg + 1) * M_PI / 4, u = sin(th);
      double g0 = u * Physical::alpha - 4 * Physical::c / u * exp(-2 * Physical::c * Rc / u) *
                  (Rc * u * Physical::alpha / 2 + pow(u * Physical::alpha, 2) / 4);
      V2 += M_PI / 4 * pow(cos(th), 2) * (1 + u * u / 2) * g0;
    }
    V2 *= K;
    UehlingSpherePotential ug2(82, Rc, 2, false, true);
    REQUIRE(V2 < 0);
    REQUIRE(ug2.V(0) == Approx(V2).epsilon(1e-12));
    REQUIRE(UehlingSpherePotential(82, Rc, 2, true).V(0) == 0);

    // At ordinary node counts: the centre agrees with the reference
    // (data/uehling_sphere_ref.py, r = 0) and joins the values just off centre
    vector<pair<int, double>> centre = {{82, -4421.0316696091921}, {14, -1581.2938391256568}};
    for (const pair<int, double> &cv : centre) {
      double Rn = (cv.first == 82 ? 7.1 : 4.0) * Physical::fm;
      UehlingSpherePotential ug(cv.first, Rn, 400, false, true);
      REQUIRE(fabs(ug.V(0) / cv.second - 1) < 1e-10);
      REQUIRE(fabs(ug.V(1e-4 * Rn) / ug.V(0) - 1) < 1e-7);
      // Below a non-zero low cutoff the same central value is returned
      double v0 = ug.V(0);
      ug.set_exp_cutoffs(1e-3, INFINITY);
      REQUIRE(ug.V(1e-12) == v0);
    }
  }

  // Kernel selection: gauss takes precedence; a point nucleus keeps its own kernel
  double R = 7.1 * Physical::fm;
  REQUIRE(UehlingSpherePotential(82, R, 100).getKernelMode() == 0);
  REQUIRE(UehlingSpherePotential(82, R, 100, true).getKernelMode() == 1);
  REQUIRE(UehlingSpherePotential(82, R, 100, false, true).getKernelMode() == 2);
  REQUIRE(UehlingSpherePotential(82, R, 100, true, true).getKernelMode() == 2);
  REQUIRE(UehlingSpherePotential(82, -1, 100, true, true).getKernelMode() == 0);

  // The options reach the potential through Atom::setUehling
  DiracAtom da = DiracAtom(82, Physical::m_mu, 208, SPHERE);
  da.setUehling(true, 100, 0, INFINITY, false, true);
  REQUIRE(da.getPotentialUehling().getKernelMode() == 2);
  da.setUehling(true, 100, 0, INFINITY, true, false);
  REQUIRE(da.getPotentialUehling().getKernelMode() == 1);
  da.setUehling(true, 100);
  REQUIRE(da.getPotentialUehling().getKernelMode() == 0);
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