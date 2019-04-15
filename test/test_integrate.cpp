#include <vector>
#include <iostream>
#include <math.h>
#include "../src/utils.hpp"
#include "../src/integrate.hpp"
#include "../src/boundary.hpp"

#include "catch/catch.hpp"

#define ERRTOL 1e-5

using namespace std;

// Integration tests
vector<double> applyFunc(double (*f)(double), vector<double> x)
{
    vector<double> y(x.size());

    for (int i = 0; i < x.size(); ++i)
    {
        y[i] = f(x[i]);
    }

    return y;
}

double trapzIntTest(double (*f)(double), double x0 = 0, double x1 = 1, int n = 200)
{
    vector<double> x, y;

    x = linGrid(x0, x1, n);
    y = applyFunc(f, x);

    return trapzInt(x, y);
}

double shootTest(double (*fQ)(double), double (*fA)(double), double (*fB)(double),
                 double x0 = 0, double x1 = 1, int n = 200, char type = 'Q', char dir = 'f')
{
    double err = 0.0;
    vector<double> x, A(n), B(n), Q(n);

    x = linGrid(x0, x1, n);
    A = applyFunc(fA, x);
    B = applyFunc(fB, x);

    // Boundary
    if (dir == 'f')
    {
        Q[0] = fQ(x[0]);
        Q[1] = fQ(x[1]);
    }
    else
    {
        Q[n - 2] = fQ(x[n - 2]);
        Q[n - 1] = fQ(x[n - 1]);
    }

    switch (type)
    {
    case 'Q':
        shootQ(Q, A, B, x[1] - x[0], dir == 'f' ? -1 : 0, dir);
        break;
    case 'N':
        shootNumerov(Q, A, B, x[1] - x[0], dir == 'f' ? -1 : 0, dir);
        break;
    default:
        break;
    }

    for (int i = 0; i < n; ++i)
    {
        err += abs(Q[i] - fQ(x[i]));
    }
    err /= n;

    return err;
}

double shootQPTest(double (*fQ)(double), double (*fP)(double),
                   double (*fAA)(double), double (*fAB)(double),
                   double (*fBA)(double), double (*fBB)(double),
                   double x0 = 0, double x1 = 1, int n = 200, char dir = 'f')
{
    double err = 0.0;
    vector<double> x, AA(n), AB(n), BA(n), BB(n), Q(n), P(n);

    x = linGrid(x0, x1, n);
    AA = applyFunc(fAA, x);
    AB = applyFunc(fAB, x);
    BA = applyFunc(fBA, x);
    BB = applyFunc(fBB, x);

    // Boundary
    if (dir == 'f')
    {
        Q[0] = fQ(x[0]);
        Q[1] = fQ(x[1]);
        P[0] = fP(x[0]);
        P[1] = fP(x[1]);
    }
    else
    {
        Q[n - 2] = fQ(x[n - 2]);
        Q[n - 1] = fQ(x[n - 1]);
        P[n - 2] = fP(x[n - 2]);
        P[n - 1] = fP(x[n - 1]);
    }

    shootQP(Q, P, AA, AB, BA, BB, x[1] - x[0], dir == 'f' ? -1 : 0, dir);

    for (int i = 0; i < n; ++i)
    {
        err += abs(Q[i] - fQ(x[i]));
        err += abs(P[i] - fP(x[i]));
    }
    err /= 2 * n;

    return err;
}

double shootPotentialTest(double (*fV)(double), double (*frho)(double),
                          double x0 = 0, double x1 = 1, int n = 200)
{
    double err = 0.0;
    vector<vector<double>> grid = logGrid(x0, x1, n);
    vector<double> rho(n), V(n);

    rho = applyFunc(frho, grid[1]);

    shootPotentialLog(V, rho, grid[0][1] - grid[0][0]);

    for (int i = 0; i < n; ++i)
    {
        err += abs(V[i] - fV(grid[1][i]));
    }
    err /= n;

    return err;
}

TEST_CASE("Trapezoidal integration", "[trapzInt]")
{
    REQUIRE(trapzIntTest(exp) == Approx(exp(1.0) - 1.0));
    REQUIRE(trapzIntTest(sin) == Approx(1.0 - cos(1.0)));
    REQUIRE(trapzIntTest(sqrt, 1.0, 4.0) == Approx(14.0 / 3.0));
    REQUIRE(trapzIntTest(cbrt, 1.0, 8.0) == Approx(45.0 / 4.0));
}

TEST_CASE("Shooting integration", "[shootQ]")
{
    /* 
        Q = tan(x)
        Q' = tan(x)^2 + 1
    */
    auto one = [](double x) { return 1.0; };
    REQUIRE(shootTest(tan, tan, one, 0, 1, 1000) < ERRTOL);
    REQUIRE(shootTest(tan, tan, one, 0, 1, 1000, 'Q', 'b') < ERRTOL);
    /* 
        Q = exp(-x**2) + x
        Q' = -2x(Q-x) + 1 = (-2x)Q+(2x^2+1)
    */
    auto fQ = [](double x) { return exp(-x * x) + x; };
    auto fA = [](double x) { return -2 * x; };
    auto fB = [](double x) { return 2 * x * x + 1; };
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000) < ERRTOL);
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000, 'Q', 'b') < ERRTOL);
}

TEST_CASE("Numerov integration", "[shootNumerov]")
{
    auto fQ = [](double x) { return sin(3 * x) + x * x; };
    auto fA = [](double x) { return -9.0; };
    auto fB = [](double x) { return 9.0 * x * x + 2; };
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000, 'N') < ERRTOL);
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000, 'N', 'b') < ERRTOL);
}

TEST_CASE("Coupled integration", "[shootQP]")
{
    auto fQ = [](double x) { return sin(x); };
    auto fP = [](double x) { return cos(x); };
    auto one = [](double x) { return 1.0; };
    auto negone = [](double x) { return -1.0; };
    auto zero = [](double x) { return 0.0; };
    REQUIRE(shootQPTest(fQ, fP, zero, one, negone, zero, 0, 1, 1000) < ERRTOL);
    REQUIRE(shootQPTest(fQ, fP, zero, one, negone, zero, 0, 1, 1000, 'b') < ERRTOL);
}

TEST_CASE("Potential integration", "[shootPotentialLog]")
{
    // Test with constant density
    auto frho = [](double x) { return 4 * M_PI * x * x; };
    auto fV = [](double x) { return 2.0 / 3.0 * M_PI * x * x; };
    REQUIRE(shootPotentialTest(fV, frho, 1e-2, 1e-2 * M_E, 1000) < ERRTOL);
    // Test with finite sphere
    auto frho_f = [](double x) { return x < 1 ? 4 * M_PI * x * x : 0; };
    auto fV_f = [](double x) { return x < 1 ? 2.0 / 3.0 * M_PI * x * x : 2.0 / 3.0 * M_PI * (3 - 2 / x); };
    REQUIRE(shootPotentialTest(fV_f, frho_f, 1e-1, 1e1, 2000) < ERRTOL);
}