#include <vector>
#include <iostream>
#include <math.h>
#include "../lib/potential.hpp"

#include "../vendor/catch/catch.hpp"

#define ERRTOL_HIGH 1e-5
#define ERRTOL_LOW 1e-3

TEST_CASE("Coulomb potential", "[CoulombPotential]")
{
    CoulombPotential cpot = CoulombPotential(1, 1.0);
    REQUIRE(cpot.V(1) == -1);
    REQUIRE(cpot.V(0) == Approx(-1.5));
    REQUIRE(cpot.V(2) == Approx(-0.5));
    REQUIRE_THROWS_AS(cpot.V(-1), invalid_argument);
}
