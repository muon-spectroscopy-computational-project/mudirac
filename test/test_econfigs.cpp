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
    REQUIRE(econf.maxn() == 3);

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

TEST_CASE("Electronic density", "[atomEdens]")
{
    ElectronicConfiguration econf("He");

    REQUIRE(econf.hydrogenicChargeDensity(1, 2) == Approx(64*exp(-4)));

    econf = ElectronicConfiguration("Li");
    
    REQUIRE(econf.hydrogenicChargeDensity(1/3.0, 3, 1.0, true) == Approx(216/9.0*exp(-2)+1.0/(9*18)*exp(-1/3.0)));
}