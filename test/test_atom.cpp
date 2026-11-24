#include "../lib/atom.hpp"
#include "../lib/hydrogenic.hpp"
#include "../lib/utils.hpp"
#include "../vendor/aixlog/aixlog.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <vector>

#include "../vendor/catch/catch.hpp"
#include "macros.hpp"

using namespace std;

TEST_CASE("Transition Matrix", "[TransitionMatrix]") {

  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::info,
                                      AixLog::Type::normal);

  TransitionMatrix tmat(-1, 2);

  REQUIRE(tmat.m1.size() == 2);
  REQUIRE(tmat.m2.size() == 4);
  REQUIRE(tmat.T.size() == tmat.m1.size());
  REQUIRE(tmat.T[0].size() == tmat.m2.size());
  REQUIRE(tmat.T[0][0] == 0);
}

TEST_CASE("Dirac Atom - basics", "[DiracAtom]") {
  double Z = 1;
  double m = 1;
  double A = 1;
  DiracAtom da = DiracAtom(Z, m, A);

  REQUIRE(da.getZ() == Z);
  REQUIRE(da.getm() == m);
  REQUIRE(da.getA() == A);
  REQUIRE(da.getmu() == Approx(effectiveMass(m, A * Physical::amu)));
  REQUIRE(da.getRestE() == da.getmu() * pow(Physical::c, 2));
}

TEST_CASE("Dirac Atom - full mass", "[DiracAtom]") {
  double Z = 1;
  double m = 1;
  double A = 1;
  DiracAtom da = DiracAtom(Z, m, A, POINT, -1, 1, 0.005, -1, false);
  REQUIRE (da.getmu() == m);
}

TEST_CASE("Dirac Atom - nuclear mass", "[DiracAtom]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::warning,
                                      AixLog::Type::all);
  // B_e(Z) = 14.4381 Z^2.39 + 1.55468e-6 Z^5.35 eV, evaluated independently
  double E_H = 27.211386245981;
  REQUIRE(Atom::electronBindingEnergy(14) * E_H == Approx(7922.668807420379).epsilon(1e-13));
  REQUIRE(Atom::electronBindingEnergy(82) * E_H == Approx(568357.7563871737).epsilon(1e-13));
  REQUIRE(Atom::electronBindingEnergy(92) * E_H == Approx(762669.862342978).epsilon(1e-13));

  // Reduced mass with M_atom - Z m_e + B_e/c^2
  int Z = 14, A = 28;
  double M_atom = getIsotopeMass(Z, A) * Physical::amu;
  double M_nuc = M_atom - Z * Physical::m_e + Atom::electronBindingEnergy(Z) / pow(Physical::c, 2);
  DiracAtom da_atom = DiracAtom(Z, Physical::m_mu, A, FERMI2);
  DiracAtom da_nuc = DiracAtom(Z, Physical::m_mu, A, FERMI2, -1, 1.0, 0.005, -1, true, true);
  REQUIRE(da_atom.getmu() == Approx(effectiveMass(Physical::m_mu, M_atom)).epsilon(1e-15));
  REQUIRE(da_nuc.getmu() == Approx(effectiveMass(Physical::m_mu, M_nuc)).epsilon(1e-15));
  REQUIRE(da_nuc.getmu() < da_atom.getmu());
  // The atomic mass is still reported by getM
  REQUIRE(da_nuc.getM() == da_atom.getM());

  // The K binding energy of muonic Si-28 drops by about 0.56 eV (1.1 ppm)
  double dE = (da_atom.getState(1, 0, true).bindingE() - da_nuc.getState(1, 0, true).bindingE()) / Physical::eV;
  REQUIRE(dE > -0.65);
  REQUIRE(dE < -0.50);

  // Without the reduced mass, the option has no effect
  DiracAtom da_full = DiracAtom(Z, Physical::m_mu, A, FERMI2, -1, 1.0, 0.005, -1, false, true);
  REQUIRE(da_full.getmu() == Physical::m_mu);
}

TEST_CASE("Dirac Atom - grid", "[DiracAtom]") {
  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::info,
                                      AixLog::Type::normal);
  double Z = 92;
  double m = 1;
  int n = 1;
  int k = -1;
  DiracAtom da = DiracAtom(Z, m);
  double E0 = hydrogenicDiracEnergy(Z, m, n, k);
  double B0 = E0 - da.getRestE();

  DiracState ds = da.initState(E0, k);

  REQUIRE(ds.grid[0] < -Z / B0);
  REQUIRE(ds.grid[ds.grid.size() - 1] > -Z / B0);

  // Must throw runtime_error for invalid out_eps
  da.out_eps = 2;
  REQUIRE_THROWS(da.gridLimits(E0, k));
}

TEST_CASE("Dirac Atom - energy search", "[DiracAtom]") {

  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::info,
                                      AixLog::Type::normal);

  double Z = 1;
  double m = 1;
  double Es1 = hydrogenicDiracEnergy(Z, m, 1),
         Es2 = hydrogenicDiracEnergy(Z, m, 2),
         Es3 = hydrogenicDiracEnergy(Z, m, 3);
  DiracAtom da = DiracAtom(Z, m);
  DiracState ds;
  TurningPoint tp;

  pair<double, double> limE = da.energyLimits();
  REQUIRE(limE.first == Approx(2 * Es1 - da.getRestE()));
  REQUIRE(limE.second == Approx(da.getRestE()));

  // Test finding the nodes interval
  int nodes;
  qnumPrincipal2Nodes(2, 0, nodes);
  LOG(INFO) << nodes << "\n";

  limE = da.energyLimits(nodes);
  ds.k = -1;

  da.convergeNodes(ds, tp, nodes, limE.first, limE.second);
  REQUIRE(ds.E > Es1);
  REQUIRE(ds.E < Es3);

  // Test finding the correct energy
  da.convergeE(ds, tp, limE.first, limE.second);
  REQUIRE(ds.E == Approx(Es2));

  // And test full convergence
  ds = da.convergeState(2, -1);
  REQUIRE(ds.E == Approx(Es2));
}

TEST_CASE("Dirac Atom - transitions", "[DiracAtom]") {
  // Tests are carried out with an ideal hydrogen atom
  // Exact values are taken from NIST database
  DiracIdealAtom daH = DiracIdealAtom(1, 1, 1, NuclearRadiusModel::SPHERE);
  TransitionMatrix tmat(-1, -1);

  // 2p3/2 => 1s1/2
  tmat = daH.getTransitionProbabilities(2, 1, true, 1, 0, true);
  REQUIRE(tmat.totalRate() * Physical::s == Approx(6.2648e+08).epsilon(1e-4));
  // 2p1/2 => 1s1/2
  tmat = daH.getTransitionProbabilities(2, 1, false, 1, 0, true);
  REQUIRE(tmat.totalRate() * Physical::s == Approx(6.2649e+08).epsilon(1e-4));

  // 3d5/2 => 2p3/2
  tmat = daH.getTransitionProbabilities(3, 2, true, 2, 1, true);
  REQUIRE(tmat.totalRate() * Physical::s == Approx(6.4651e+07).epsilon(1e-4));

  // 3d3/2 => 2p3/2
  tmat = daH.getTransitionProbabilities(3, 2, false, 2, 1, true);
  REQUIRE(tmat.totalRate() * Physical::s == Approx(1.0775e+07).epsilon(1e-4));

  // Iron-like atom
  DiracAtom daFe = DiracAtom(26, 1, 56, NuclearRadiusModel::SPHERE);

  // 3d3/2 => 3p1/2
  tmat = daFe.getTransitionProbabilities(3, 2, false, 3, 1, false);
  REQUIRE(
    tmat.totalRate() * Physical::s ==
    Approx(1.31e7).epsilon(
      3e-2)); // Precision is not strong here... possibly needs improvement
}