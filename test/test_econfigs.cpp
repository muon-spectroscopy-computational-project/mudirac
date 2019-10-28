#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "../lib/econfigs.hpp"

#include "../vendor/catch/catch.hpp"

using namespace std;

TEST_CASE("Atomic electronic configurations", "[atomEconfig]")
{
    ElectronicConfiguration econf("1s2");
    REQUIRE(econf.getPopulation(1, 0) == 2); // 1s2

    econf = ElectronicConfiguration("[He] 2s2");
    REQUIRE(econf.getPopulation(1, 0) == 2); // 1s2
    REQUIRE(econf.getPopulation(2, 0) == 2); // 2s2

    econf = ElectronicConfiguration("Na");
    REQUIRE(econf.getPopulation(1, 0) == 2); // 1s2
    REQUIRE(econf.getPopulation(2, 0) == 2); // 2s2
    REQUIRE(econf.getPopulation(2, 1) == 6); // 2p6
    REQUIRE(econf.getPopulation(3, 0) == 1); // 3s1

    // Exceptions to the Madelung rule
    econf = ElectronicConfiguration("Cu");
    REQUIRE(econf.getPopulation(3, 2) == 10); // 3d10

    // Failure cases
    REQUIRE_THROWS(ElectronicConfiguration("nonsense"));
    REQUIRE_THROWS(ElectronicConfiguration("1s2 1s4"));
    REQUIRE_THROWS(ElectronicConfiguration("[He] 1s2"));
    REQUIRE_THROWS(ElectronicConfiguration("[He] [Ne]"));
    REQUIRE_THROWS(ElectronicConfiguration("5a2"));
}