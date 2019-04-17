/**
 * Test functions for utils.cpp
 */

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "../src/utils.hpp"
#include "data/genlag.h"

#include "catch/catch.hpp"

bool compareLaguerre(double tol = 1e-5)
{
    // Compare Laguerre polynomials with reference table

    int val;
    int pN, xN;
    string tag, error;

    int n, alpha;
    double x, ydat, ycalc;

    pN = sizeof(laguerreParams) / (2 * sizeof(int));
    xN = sizeof(laguerreValues) / ((pN + 1) * sizeof(double));

    // Now read everything else
    for (int i = 0; i < xN; ++i)
    {
        x = laguerreValues[i][0];
        for (int j = 0; j < pN; ++j)
        {
            n = laguerreParams[j][0];
            alpha = laguerreParams[j][1];
            ydat = laguerreValues[i][j + 1];
            ycalc = genLaguerrePoly(x, n, alpha);
            if (abs(ydat - ycalc) > tol)
            {
                throw("Polynomials don't match for x = " + to_string(x) +
                      " n = " + to_string(n) +
                      " alpha = " + to_string(alpha));
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

TEST_CASE("Splitting strings", "[splitstring]")
{
    // Basic
    CHECK(splitString("a,b,c", ",") == vector<string>{"a", "b", "c"});
    // Merge
    CHECK(splitString("a  b c") == vector<string>{"a", "", "b", "c"});
    CHECK(splitString("a  b c", " ", true) == vector<string>{"a", "b", "c"});
    // Multiple separators
    CHECK(splitString("a,b c", ", ") == vector<string>{"a", "b", "c"});
    // Separators at the beginning
    CHECK(splitString("  a,b c", ", ", true) == vector<string>{"a", "b", "c"});
}