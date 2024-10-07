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

TEST_CASE("Transition Matrix", "[TransitionMatrix]")
{

  AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::info,
                                      AixLog::Type::normal);

  TransitionMatrix tmat(-1, 2);

  REQUIRE(tmat.m1.size() == 2);
  REQUIRE(tmat.m2.size() == 4);
  REQUIRE(tmat.T.size() == tmat.m1.size());
  REQUIRE(tmat.T[0].size() == tmat.m2.size());
  REQUIRE(tmat.T[0][0] == 0);
}

TEST_CASE("Dirac Atom - basics", "[DiracAtom]")
{
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

TEST_CASE("Dirac Atom - grid", "[DiracAtom]")
{
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

TEST_CASE("Dirac Atom - energy search", "[DiracAtom]")
{

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

TEST_CASE("Dirac Atom - transitions", "[DiracAtom]")
{
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