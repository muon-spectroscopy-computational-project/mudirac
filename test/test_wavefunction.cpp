/**
 * Test functions for wavefunction.cpp (operations with wavefunctions, bra-kets and the like)
 */

#include <iostream>
#include <vector>
#include "../lib/hydrogenic.hpp"
#include "../lib/utils.hpp"
#include "../lib/wavefunction.hpp"

#include "../vendor/catch/catch.hpp"

using namespace std;

TEST_CASE("Bra-ket on logarithmic grid", "[braKetLog]")
{
    // Define a grid
    double dx = 0.01;
    vector<vector<double>> grid = logGrid(1, dx, -500, 500);

    vector<double> s1 = hydrogenicSchroWavefunction(grid[1], 1, 1, 1, 0);
    vector<double> s2 = hydrogenicSchroWavefunction(grid[1], 1, 1, 2, 0);

    REQUIRE(braKetLog(s1, s1, grid[1], dx) == Approx(1));
    REQUIRE(braKetLog(s2, s2, grid[1], dx) == Approx(1));
    REQUIRE(braKetLog(s1, s2, grid[1], dx) + 1 == Approx(1));

    vector<vector<double>> s1d = hydrogenicDiracWavefunction(grid[1], 1, 1, 1, -1);
    vector<vector<double>> s2d = hydrogenicDiracWavefunction(grid[1], 1, 1, 2, -1);

    REQUIRE(braKetLog(s1d, s1d, grid[1], dx) == Approx(1));
    REQUIRE(braKetLog(s2d, s2d, grid[1], dx) == Approx(1));
    REQUIRE(braKetLog(s1d, s2d, grid[1], dx) + 1 == Approx(1));
}

TEST_CASE("Expectation values on logarithmic grid", "[braOpKetLog]")
{
    // Define a grid
    double dx = 0.01;
    vector<vector<double>> grid = logGrid(1, dx, -500, 500);

    vector<double> s1 = hydrogenicSchroWavefunction(grid[1], 1, 1, 1, 0);
    vector<double> s2 = hydrogenicSchroWavefunction(grid[1], 1, 1, 2, 0);

    REQUIRE(braOpKetLog(s1, grid[1], s1, grid[1], dx) == Approx(1.5));
}