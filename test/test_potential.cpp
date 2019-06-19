#include <vector>
#include <iostream>
#include <math.h>
#include "../lib/potential.hpp"
#include "../lib/atom.hpp"

#include "../vendor/catch/catch.hpp"

#define ERRTOL_HIGH 1e-5
#define ERRTOL_LOW 1e-3

TEST_CASE("Coulomb potential", "[CoulombSpherePotential]")
{
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

TEST_CASE("Uehling correction to Coulomb potential", "[CoulombSpherePotential]")
{
    UehlingSpherePotential cpot = UehlingSpherePotential(1, 1e-4, 100);

    // Test static functions
    REQUIRE(cpot.ukernel_r_greater(0.5, 1, 0) == Approx(0));
    REQUIRE(cpot.ukernel_r_smaller(0.5, 1, 0) == Approx(0));
    REQUIRE(cpot.ukernel_r_smaller(0.5, 0, 1) == Approx(0));
    REQUIRE(cpot.ukernel_point(0.5, 0) == Approx(2));

    // More complex cases
    REQUIRE(cpot.ukernel_r_greater(0.5, 0.5 / (2 * Physical::c), 0.5 / (2 * Physical::c)) == Approx(0.125 / pow(Physical::c * M_E, 2)));
    REQUIRE(cpot.ukernel_r_smaller(0.5, 0.5 / (2 * Physical::c), 0.5 / (2 * Physical::c)) ==
            Approx(2 * sinh(1) * 0.125 / pow(Physical::c, 2) * (0.5 - exp(-1))));
    REQUIRE(cpot.ukernel_point(0.5, 0.5 / (2 * Physical::c)) == Approx(2 * exp(-1)));

    // Now the actual potential
    cpot = UehlingSpherePotential(92, Atom::sphereNuclearModel(238), 500);
    REQUIRE(cpot.V(0.01) == Approx(-0.1796).epsilon(1e-3));
    REQUIRE(cpot.V(5e-5) == Approx(-4438.7471).epsilon(1e-3));
    // Point-like atom
    cpot = UehlingSpherePotential(92, -1, 2000);
    REQUIRE(cpot.V(0.01) == Approx(-0.1795).epsilon(1e-3));
    REQUIRE(cpot.V(5e-5) == Approx(-10225.5698).epsilon(1e-3));
}