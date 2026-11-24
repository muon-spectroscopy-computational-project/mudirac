#include "../lib/atom.hpp"
#include "../lib/config.hpp"
#include "../lib/optimisation.hpp"
#include "../lib/potential.hpp"
#include "datapath.h"
#include <fstream>
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

TEST_CASE("Uehling potential of a Fermi 2-term density", "[UehlingFermi2Potential]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::warning,
                                      AixLog::Type::all);
  double t = 2.3 * Physical::fm;

  // Reference values from data/fermi2_uehling_ref.py: the direct convolution of
  // the density with the point-charge Uehling potential, by mpmath (an
  // independent route from the superposition of uniform spheres used here).
  // c from MuDirac's rule for R = 7.102018661385413 fm (Pb-208) and
  // R = 4.031001066732679 fm (Si-28), t = 2.3 fm; r = x*c
  struct UehlingRef {
    int Z;
    double c, x, V;
  };
  vector<UehlingRef> refs = {
    {82, 0.00012553558838635206, 0.1, -4578.1694452265975},
    {82, 0.00012553558838635206, 0.5, -4154.1161227326006},
    {82, 0.00012553558838635206, 1.0, -2803.6033898240004},
    {82, 0.00012553558838635206, 1.5, -1594.7363873226436},
    {82, 0.00012553558838635206, 3.0, -567.88285579431864},
    {82, 0.00012553558838635206, 10.0, -72.436314577765446},
    {14, 5.958037111956458e-05, 0.1, -1767.4186414571072},
    {14, 5.958037111956458e-05, 0.5, -1615.5648472980248},
    {14, 5.958037111956458e-05, 1.0, -1178.0194080392392},
    {14, 5.958037111956458e-05, 1.5, -748.08524136685443},
    {14, 5.958037111956458e-05, 3.0, -288.03585595235353},
    {14, 5.958037111956458e-05, 10.0, -46.550884589768964}
  };
  for (const UehlingRef &ref : refs) {
    UehlingFermi2Potential vu(ref.Z, ref.c, t, 100, 12);
    UehlingFermi2Potential vu_direct(ref.Z, ref.c, t, 100, 12, true);
    double r = ref.x * ref.c;
    double v = vu.V(r);
    REQUIRE(fabs(v / ref.V - 1) < 1e-11);
    // The precomputed panel sums reproduce the plain double sum
    REQUIRE(fabs(v / vu_direct.V(r) - 1) < 2e-12);
    // Cached values are returned as computed
    REQUIRE(vu.V(r) == v);
    REQUIRE(vu.cacheSize() == 1);
  }

  // h(x) = (x - 1) + (x + 1) exp(-2x): series below 1, closed form above
  for (double x : {0.3, 0.6, 0.9}) {
    double h = (x - 1) + (x + 1) * exp(-2 * x);
    REQUIRE(UehlingFermi2Potential::hOut(x) == Approx(h).epsilon(1e-13));
  }
  REQUIRE(UehlingFermi2Potential::hOut(1 - 1e-12) == Approx(UehlingFermi2Potential::hOut(1)).epsilon(1e-10));
  // Leading term 2x^3/3 for small x, where the closed form cancels
  REQUIRE(UehlingFermi2Potential::hOut(1e-4) == Approx(2e-12 / 3 * (1 - 1e-4)).epsilon(1e-7));

  REQUIRE_THROWS_AS(UehlingFermi2Potential(82, 1e-4, t, 0, 12), invalid_argument);
  REQUIRE_THROWS_AS(UehlingFermi2Potential(82, 1e-4, t, 100, 0), invalid_argument);
  REQUIRE_THROWS_AS(UehlingFermi2Potential(82, 1e-4, -t, 100, 12), invalid_argument);
  REQUIRE_THROWS_AS(UehlingFermi2Potential(82, 1e-4, t).V(0.0), invalid_argument);
}

TEST_CASE("Uehling potential of a Fermi 2-term density in an atom", "[DiracAtom]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::warning,
                                      AixLog::Type::all);

  // The Uehling potential follows the Coulomb potential's c and t
  DiracAtom da = DiracAtom(82, Physical::m_mu, 208, FERMI2, 7.102018661385413);
  da.setUehlingDensity(true);
  da.setUehling(true, 100);
  CoulombFermi2Potential *vc = (CoulombFermi2Potential *)da.getPotentialCoulomb();
  REQUIRE(da.getUehlingDensityFermi2());
  REQUIRE(da.getPotentialUehlingFermi2().getc() == vc->getc());
  REQUIRE(da.getPotentialUehlingFermi2().getT() == vc->getT());
  da.setFermi2(2.0 * Physical::fm);
  vc = (CoulombFermi2Potential *)da.getPotentialCoulomb();
  REQUIRE(da.getPotentialUehlingFermi2().getc() == vc->getc());
  REQUIRE(da.getPotentialUehlingFermi2().getT() == 2.0 * Physical::fm);
  da.setFermi2Femto(6.6, 2.3, CT);
  vc = (CoulombFermi2Potential *)da.getPotentialCoulomb();
  REQUIRE(da.getPotentialUehlingFermi2().getc() == vc->getc());
  REQUIRE(vc->getc() == 6.6 * Physical::fm);
  // ...and is the one Atom::getV adds
  double r = 5e-5;
  UehlingFermi2Potential vu(82, vc->getc(), vc->getT(), 100, 12);
  REQUIRE(da.getV(r) == vc->V(r) + vu.V(r));

  // Refusals: a nucleus that is not FERMI2, no shell nodes, a low cutoff
  DiracAtom da_s = DiracAtom(82, Physical::m_mu, 208, SPHERE);
  // (each refusal is identified by its message, so that no guard hides behind another)
  REQUIRE_THROWS_WITH(da_s.setUehlingDensity(true),
                      "The Uehling potential of the Fermi 2-term density needs a FERMI2 nuclear model");
  REQUIRE_THROWS_WITH(da.setUehlingDensity(true, 0),
                      "The Uehling potential of the Fermi 2-term density needs at least 1 shell node");
  REQUIRE_THROWS_WITH(da.setUehling(true, 100, 1e-3),
                      "The Uehling potential of the Fermi 2-term density does not support a low cutoff");
  // A FERMI2 atom without a nuclear radius (no isotope) has no density to use
  DiracAtom da_noR = DiracAtom(82, Physical::m_mu, -1, FERMI2);
  da_noR.setUehlingDensity(true);
  REQUIRE_THROWS_WITH(da_noR.setUehling(true, 100),
                      "The Uehling potential of the Fermi 2-term density needs a finite FERMI2 nucleus");
}

// Write an input file into the test data folder and build its atom
static DiracAtom makeAtomFromText(string name, string text) {
  string path = string(CURRENT_DATAPATH) + "/data/" + name;
  ofstream out(path);
  out << text;
  out.close();
  MuDiracInputFile config;
  config.parseFile(path);
  return config.makeAtom();
}

TEST_CASE("Uehling potential of a Fermi 2-term density from the input file", "[MuDiracInputFile]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::fatal,
                                      AixLog::Type::all);
  string base = "element: Pb\nisotope: 208\nuehling_correction: T\n";

  DiracAtom da = makeAtomFromText("uehling_fermi2_ok.in", base + "nuclear_model: FERMI2\nuehling_density: fermi2\nuehling_shell_nodes: 16\n");
  REQUIRE(da.getUehlingDensityFermi2());
  REQUIRE(da.getPotentialUehlingFermi2().getSnodes() == 16);
  da = makeAtomFromText("uehling_sphere_ok.in", base + "nuclear_model: FERMI2\n");
  REQUIRE(!da.getUehlingDensityFermi2());

  REQUIRE_THROWS_WITH(makeAtomFromText("uehling_bad_density.in", base + "nuclear_model: FERMI2\nuehling_density: gauss\n"),
                      "Invalid uehling_density parameter in input file");
  REQUIRE_THROWS_WITH(makeAtomFromText("uehling_bad_model.in", base + "nuclear_model: SPHERE\nuehling_density: FERMI2\n"),
                      "uehling_density FERMI2 needs nuclear_model FERMI2");
  REQUIRE_THROWS_WITH(makeAtomFromText("uehling_bad_lowcut.in", base + "nuclear_model: FERMI2\nuehling_density: FERMI2\nuehling_lowcut: 1e-3\n"),
                      "uehling_density FERMI2 does not support uehling_lowcut");
  REQUIRE_THROWS_WITH(makeAtomFromText("uehling_bad_nodes.in", base + "nuclear_model: FERMI2\nuehling_density: FERMI2\nuehling_shell_nodes: 0\n"),
                      "The Uehling potential of the Fermi 2-term density needs at least 1 shell node");
}

TEST_CASE("Uehling potential of a Fermi 2-term density, refused settings", "[DiracAtom]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::fatal,
                                      AixLog::Type::all);
  double r = 5e-5;
  string lowcut_msg = "The Uehling potential of the Fermi 2-term density does not support a low cutoff";

  // A refused setUehling leaves the atom as it was (here: no Uehling term)
  DiracAtom da1 = DiracAtom(82, Physical::m_mu, 208, FERMI2, 7.102018661385413);
  da1.setUehlingDensity(true);
  REQUIRE_THROWS_WITH(da1.setUehling(true, 100, 1e-3), lowcut_msg);
  REQUIRE(!da1.getUehling());
  REQUIRE(std::isfinite(da1.getV(r)));
  REQUIRE(da1.getV(r) == da1.getPotentialCoulomb()->V(r));

  // A refused setUehlingDensity leaves the uniform-sphere Uehling term in place
  DiracAtom da2 = DiracAtom(82, Physical::m_mu, 208, FERMI2, 7.102018661385413);
  da2.setUehling(true, 100, 1e-3);
  REQUIRE_THROWS_WITH(da2.setUehlingDensity(true), lowcut_msg);
  REQUIRE(!da2.getUehlingDensityFermi2());
  UehlingSpherePotential vs(82, 7.102018661385413 * Physical::fm, 100);
  vs.set_exp_cutoffs(1e-3, INFINITY);
  REQUIRE(std::isfinite(da2.getV(r)));
  REQUIRE(da2.getV(r) == da2.getPotentialCoulomb()->V(r) + vs.V(r));

  // A refused rebuild (negative thickness) keeps the previous Coulomb and
  // Uehling potentials together, through both Fermi 2-term setters
  DiracAtom da3 = DiracAtom(82, Physical::m_mu, 208, FERMI2, 7.102018661385413);
  da3.setUehlingDensity(true);
  da3.setUehling(true, 100);
  da3.setFermi2(2.3 * Physical::fm);
  double V3 = da3.getV(r);
  double c3 = ((CoulombFermi2Potential *)da3.getPotentialCoulomb())->getc();
  REQUIRE_THROWS_AS(da3.setFermi2(-0.5 * Physical::fm), invalid_argument);
  REQUIRE(da3.getV(r) == V3);
  REQUIRE(((CoulombFermi2Potential *)da3.getPotentialCoulomb())->getc() == c3);
  REQUIRE(da3.fermi2.t == 2.3);
  REQUIRE_THROWS_AS(da3.setFermi2Femto(6.6, -0.5, CT), invalid_argument);
  REQUIRE(da3.getV(r) == V3);
  REQUIRE(da3.fermi2.t == 2.3);
}

TEST_CASE("Uehling potential of a Fermi 2-term density, zero thickness", "[UehlingFermi2Potential]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::fatal,
                                      AixLog::Type::all);
  // t = 0 is the uniform sphere of radius c. Reference values from
  // data/fermi2_uehling_ref.py (the same convolution for a uniform sphere)
  double c = 6.6 * Physical::fm;
  UehlingFermi2Potential v0(82, c, 0, 100, 12);
  vector<pair<double, double>> sphere_ref = {
    {0.01, -4864.2394117990238}, {0.5, -4393.6728358296788},
    {1.5, -1603.2485490344265}, {4.0, -364.73278871879436}
  };
  for (const pair<double, double> &sr : sphere_ref) {
    double r = sr.first * c;
    double v = v0.V(r);
    REQUIRE(std::isfinite(v));
    REQUIRE(fabs(v / sr.second - 1) < 1e-10);
    // ...and the limit of t -> 0+
    double d3 = UehlingFermi2Potential(82, c, 1e-3 * Physical::fm, 100, 12).V(r) / v - 1;
    double d4 = UehlingFermi2Potential(82, c, 1e-4 * Physical::fm, 100, 12).V(r) / v - 1;
    REQUIRE(fabs(d3) < 1e-7);
    REQUIRE(fabs(d4) < fabs(d3));
  }
  REQUIRE_THROWS_AS(UehlingFermi2Potential(82, c, -1e-3 * Physical::fm, 100, 12), invalid_argument);

  // The fit's cost function at the zero-thickness boundary, in both coordinate
  // systems (CT: t = 0; POLAR: theta = 0), gives a finite residual
  for (Fermi2CoordinateSystem cs : {CT, POLAR}) {
    DiracAtom da = DiracAtom(14, Physical::m_mu, 28, FERMI2);
    da.setUehlingDensity(true);
    da.setUehling(true, 100);
    da.setFermi2(2.3 * Physical::fm);
    da.transqnums = {TransLineSpec{1, 2, 0, 1, true, false}};
    da.setExpOptData(cs, {"K1-L2"}, {400000.0}, {1.0});
    CostFunctor cost(da);
    double c1 = (cs == CT) ? 3.0 : 3.1, c2 = 0, residual = NAN;
    REQUIRE(cost(&c1, &c2, &residual));
    REQUIRE(std::isfinite(residual));
  }
}
