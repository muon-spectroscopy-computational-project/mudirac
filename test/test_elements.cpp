#include <stdio.h>
#include <iostream>
#include <fstream>
#include "../lib/elements.hpp"

#include "catch/catch.hpp"

using namespace std;

TEST_CASE("Atomic data", "[atomData]")
{
    REQUIRE(getElementZ("C") == 6);
    REQUIRE(getElementSymbol(8) == "O");
}