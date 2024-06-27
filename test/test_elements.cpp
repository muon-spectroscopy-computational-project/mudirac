#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "../lib/elements.hpp"

#include "../vendor/catch/catch.hpp"

using namespace std;

TEST_CASE("Atomic data", "[atomData]")
{
    REQUIRE(getElementZ("C") == 6);
    REQUIRE(getElementSymbol(8) == "O");
    REQUIRE(getElementMainIsotope(6) == 12);
    REQUIRE(getIsotopeMass("C") == Approx(12.0));
    REQUIRE(getIsotopeMass(6) == Approx(12.0));
    REQUIRE(getIsotopeSpin(1, 2) == Approx(1.0));
    REQUIRE_THROWS_AS(getElementZ("Kk"), invalid_argument);
    REQUIRE_THROWS_AS(getIsotopeMass(240), invalid_argument);
    CHECK(getAllIsotopes("H") == vector<int>{1, 2, 3, 4, 5, 6, 7});
}