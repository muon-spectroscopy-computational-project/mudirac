/**
 * Test functions for utils.cpp
 */

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "../lib/utils.hpp"
#include "data/genlag.h"

#include "../vendor/catch/catch.hpp"
#include "macros.hpp"

void compareLaguerre(double tol = 1e-5)
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

TEST_CASE("Vector contains", "[vectorContains]")
{
    vector<string> v = {"H", "He", "Li", "Be"};

    CHECK(vectorContains<string>(v, "Li"));
    CHECK(!vectorContains<string>(v, "Ne"));
}

TEST_CASE("Vector operations", "[vectorOperation]")
{
    vector<double> v1(3, 2), v2(3, 3);

    CHECK(vectorOperation(v1, v2, '*') == vector<double>{6, 6, 6});
    CHECK(vectorOperation(v1, v2, '+') == vector<double>{5, 5, 5});

    CHECK(vectorOperation(v1, 3.0, '^') == vector<double>{8, 8, 8});
}

TEST_CASE("Logarithmic grids", "[logGrid]")
{
    // Test center + step form
    vector<vector<double>> lgrid_01 = logGrid(2.5 / 32.0, 2.5 * 32.0, 11);
    vector<vector<double>> lgrid_cs = logGrid(2.5, log(2.0), -5, 5);

    CHECK_APPROX(lgrid_01[1], lgrid_cs[1], 1e-5);
}

TEST_CASE("Sinc", "[sinc]")
{
    REQUIRE(sinc(0) == 1);
    REQUIRE(1 - sinc(M_PI) == Approx(1.0));
    REQUIRE(sinc(M_PI / 2.0) == Approx(2.0 / M_PI));
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
    CHECK(splitString("abc") == vector<string>{"abc"});
    CHECK(splitString("a,b,c", ",") == vector<string>{"a", "b", "c"});
    // Merge
    CHECK(splitString("a  b c") == vector<string>{"a", "", "b", "c"});
    CHECK(splitString("a  b c", " ", true) == vector<string>{"a", "b", "c"});
    // Multiple separators
    CHECK(splitString("a,b c", ", ") == vector<string>{"a", "b", "c"});
    // Separators at the beginning
    CHECK(splitString("  a,b c", ", ", true) == vector<string>{"a", "b", "c"});
    // Maximum count
    CHECK(splitString("a,b,c", ",", false, 1) == vector<string>{"a", "b,c"});
}

TEST_CASE("Stripping strings", "[stripstring]")
{
    REQUIRE(stripString("   abc  ") == "abc");
    REQUIRE(stripString("abc") == "abc");
    REQUIRE(stripString("   \t").size() == 0);
}

TEST_CASE("Uppercase strings", "[upperstring]")
{
    REQUIRE(upperString("yes") == "YES");
}

TEST_CASE("Analysing IUPAC state labels", "[parseIupacState]")
{
    int n, l;
    bool s;

    parseIupacState("K", n, l, s);
    REQUIRE(n == 1);
    REQUIRE(l == 0);
    parseIupacState("L2", n, l, s);
    REQUIRE(n == 2);
    REQUIRE(l == 1);
    REQUIRE(s == false);
    parseIupacState("M5", n, l, s);
    REQUIRE(n == 3);
    REQUIRE(l == 2);
    REQUIRE(s == true);

    vector<int> nrange, lrange;
    vector<bool> srange;

    parseIupacRange("K1:M2", nrange, lrange, srange);
    CHECK(nrange == vector<int>{1, 2, 2, 2, 3, 3});
    CHECK(lrange == vector<int>{0, 0, 1, 1, 0, 1});
    CHECK(srange == vector<bool>{false, false, false, true, false, false});
    parseIupacRange("L3", nrange, lrange, srange);
    CHECK(nrange == vector<int>{2});
    CHECK(lrange == vector<int>{1});
    CHECK(srange == vector<bool>{true});

    // String printing
    CHECK(printIupacState(1, 0, false) == "K1");
    CHECK(printIupacState(1, 0, true) == "K1");
    CHECK(printIupacState(2, 1, false) == "L2");
    CHECK(printIupacState(3, 2, true) == "M5");

    // Exceptions
    REQUIRE_THROWS(parseIupacState("B2", n, l, s));
    REQUIRE_THROWS(parseIupacState("L4", n, l, s));
    REQUIRE_THROWS(parseIupacState("M", n, l, s));
    REQUIRE_THROWS(parseIupacRange("K1:L2:M3", nrange, lrange, srange));
    REQUIRE_THROWS(printIupacState(-1, 0, false));
    REQUIRE_THROWS(printIupacState(1, 2, false));
}

TEST_CASE("Quantum numbers and Clebsch-Gordan coefficients", "[quantumNumbersAndCG]")
{
    int n, l, m, k;
    bool s;

    k = -1;
    qnumDirac2Schro(k, l, s);
    REQUIRE(((l == 0) && s));

    k = 2;
    qnumDirac2Schro(k, l, s);
    REQUIRE(((l == 2) && !s));

    l = 1;
    s = true;
    qnumSchro2Dirac(l, s, k);
    REQUIRE(k == -2);

    // Clebsch-Gordan coefficients
    REQUIRE(cgCoeff(-1, 0.5, true) == Approx(1));
    REQUIRE(cgCoeff(-1, -0.5, true) == Approx(0));
    REQUIRE(cgCoeff(-3, 2.5, true) == Approx(1));
    REQUIRE(cgCoeff(-3, 1.5, false) == Approx(sqrt(0.2)));
    REQUIRE(cgCoeff(1, 0.5, true) == Approx(-sqrt(1.0 / 3.0)));
    REQUIRE(cgCoeff(2, 0.5, false) == Approx(sqrt(3.0 / 5.0)));
    REQUIRE(cgCoeff(2, 0.5, true) == Approx(-sqrt(2.0 / 5.0)));
}