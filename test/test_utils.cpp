/**
 * Test functions for utils.cpp
 */

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "../src/utils.hpp"

#include "catch/catch.hpp"

bool compareLaguerre(double tol = 1e-5)
{
    // Compare Laguerre polynomials with reference table

    int val;
    string tag, error;
    vector<int> n, alpha;

    double x, ydat, ycalc;

    ifstream datFile;

    datFile.open("test/data/genlag.dat");

    if (!datFile)
    {
        throw "Reference table not found.";
    }

    // Read the header

    while (datFile >> tag && tag.compare("##"))
    {
        datFile >> val;
        n.push_back(val);
        datFile >> val;
        alpha.push_back(val);
    }
    datFile.ignore(100, '\n');

    if (n.size() == 0)
    {
        throw "No parameters found in reference table.";
    }

    // Now read everything else
    while (datFile >> x)
    {
        for (int i = 0; i < n.size(); ++i)
        {
            datFile >> ydat;
            ycalc = genLaguerrePoly(x, n[i], alpha[i]);
            if (abs(ydat - ycalc) > tol)
            {
                throw("Polynomials don't match for x = " + to_string(x) +
                      " n = " + to_string(n[i]) +
                      " alpha = " + to_string(alpha[i]));
            }
        }
    }
}

bool countTest(int N)
{
    // Generate a function with N nodes
    vector<double> f;

    for (int i = 0; i < 1000; ++i)
    {
        f.push_back(cos(N * M_PI * 1e-3 * i));
    }

    return (N == countNodes(f));
}

TEST_CASE("Factorials", "[factorial]")
{
    REQUIRE(factorial(0) == 1);
    REQUIRE(factorial(1) == 1);
    REQUIRE(factorial(4) == 24);
    REQUIRE(factorial(6) == 720);
    REQUIRE_THROWS(factorial(-1));
}

TEST_CASE("Generalised Laguerre polynomials", "[genlagpoly]")
{
    REQUIRE_NOTHROW(compareLaguerre());
}

TEST_CASE("Node counting", "[countnodes]")
{
    REQUIRE(countTest(1));
    REQUIRE(countTest(2));
    REQUIRE(countTest(5));
    REQUIRE(countTest(10));
}