#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <tuple>
#include <limits>
#include "../lib/atom.hpp"
#include "../lib/hydrogenic.hpp"
#include "../vendor/aixlog/aixlog.hpp"

#include "../vendor/catch/catch.hpp"

using namespace std;

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