/**
 * Test functions for utils.cpp
 */

#include <iostream>
#include "../src/utils.hpp"

#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"

TEST_CASE( "Factorials", "[factorial]") {
    REQUIRE( factorial(0) == 1);
    REQUIRE( factorial(1) == 1);
    REQUIRE( factorial(4) == 24);
    REQUIRE_THROWS( factorial(-1));
}