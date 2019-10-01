#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <tuple>
#include <limits>
#include "../lib/atom.hpp"
#include "../lib/utils.hpp"
#include "../lib/hydrogenic.hpp"
#include "../vendor/aixlog/aixlog.hpp"

#include "../vendor/catch/catch.hpp"
#include "macros.hpp"

using namespace std;

TEST_CASE("Transition Matrix", "[TransitionMatrix]")
{
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
    AixLog::Log::init<AixLog::SinkCout>(AixLog::Severity::trace, AixLog::Type::normal);

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
    double minE = limE.first;
    double maxE = limE.second;
    qnumPrincipal2Nodes(2, 0, nodes);

    ds.k = -1;
    da.convergeNodes(ds, tp, nodes, minE, maxE);
    REQUIRE(ds.E > Es1);
    REQUIRE(ds.E < Es3);

    // Test finding the correct energy
    da.convergeE(ds, tp, minE, maxE);
    REQUIRE(ds.E == Approx(Es2));

    // And test full convergence
    ds = da.convergeState(2, -1);
    REQUIRE(ds.E == Approx(Es2));
}

TEST_CASE("Dirac Atom - transitions", "[DiracAtom]")
{
    double Z = 1;
    double m = 1;

    DiracAtom da = DiracAtom(Z, m);

    TransitionMatrix tmat = da.getTransitionProbabilities(2, 1, true, 1, 0, true);

    cout << '\t';
    for (int m2 = 0; m2 < tmat.m2.size(); ++m2) {
        cout << tmat.m2[m2] << '\t'<< '\t';
    }
    cout << '\n';
    for (int m1 = 0; m1 < tmat.m1.size(); ++m1) {
        cout << tmat.m1[m1] << '\t';
        for (int m2 = 0; m2 < tmat.m2.size(); ++m2) {
            cout << tmat.T[m1][m2]*Physical::s << '\t';
        }
        cout << '\n';
    }    
}