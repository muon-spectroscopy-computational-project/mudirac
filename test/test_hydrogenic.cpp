/**
 * Test functions for hydrogenic.cpp (known solutions for hydrogen-like atoms)
 */

#include <iostream>
#include <vector>
#include "../lib/hydrogenic.hpp"
#include "../lib/utils.hpp"
#include "../lib/integrate.hpp"

#include "../vendor/catch/catch.hpp"

using namespace std;

TEST_CASE("Schroedinger atom energies", "[hydrogenicSchroEnergy]")
{
    REQUIRE(hydrogenicSchroEnergy(1, 1, 1) == Approx(-0.5));
    REQUIRE(hydrogenicSchroEnergy(1, 1, 2) == Approx(-0.5 / 4));
    REQUIRE(hydrogenicSchroEnergy(1, 1, 3) == Approx(-0.5 / 9));
}

TEST_CASE("Schroedinger atom wavefunctions", "[hydrogenicSchroWavefunction]")
{
    // Define a grid
    vector<vector<double>> grid = logGrid(1, 1e-2, -500, 500);

    vector<double> s1 = hydrogenicSchroWavefunction(grid[1], 1, 1, 1, 0);
    vector<double> s2 = hydrogenicSchroWavefunction(grid[1], 1, 1, 2, 0);

    vector<double> rho(s1.size());

    for (int i = 0; i < rho.size(); ++i)
    {
        rho[i] = pow(s1[i], 2.0) * grid[1][i];
    }

    REQUIRE(trapzInt(grid[0], rho) == Approx(1));

    for (int i = 0; i < rho.size(); ++i)
    {
        rho[i] = pow(s2[i], 2.0) * grid[1][i];
    }

    REQUIRE(trapzInt(grid[0], rho) == Approx(1));

    for (int i = 0; i < rho.size(); ++i)
    {
        rho[i] = s1[i] * s2[i] * grid[1][i];
    }

    REQUIRE(abs(trapzInt(grid[0], rho)) < 1e-5);
}

TEST_CASE("Dirac atom wavefunctions", "[hydrogenicDiracWavefunction]")
{
    // Define a grid
    vector<vector<double>> grid = logGrid(1, 0.01, -500, 500);

    vector<vector<double>> s1 = hydrogenicDiracWavefunction(grid[1], 1, 1, 1, -1);
    vector<vector<double>> s2 = hydrogenicDiracWavefunction(grid[1], 1, 1, 2, -1);

    vector<double> rho(s1[0].size());

    for (int i = 0; i < rho.size(); ++i)
    {
        rho[i] = (pow(s1[0][i], 2.0) + pow(s1[1][i], 2.0)) * grid[1][i];
    }

    REQUIRE(trapzInt(grid[0], rho) == Approx(1));

    for (int i = 0; i < rho.size(); ++i)
    {
        rho[i] = (pow(s2[0][i], 2.0) + pow(s2[1][i], 2.0)) * grid[1][i];
    }

    REQUIRE(trapzInt(grid[0], rho) == Approx(1));

    for (int i = 0; i < rho.size(); ++i)
    {
        rho[i] = (s1[0][i] * s2[0][i] + s1[1][i] * s2[1][i]) * grid[1][i];
    }

    REQUIRE(abs(trapzInt(grid[0], rho)) < 1e-5);
}