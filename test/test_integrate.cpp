#include <vector>
#include <iostream>
#include <math.h>
#include "../src/utils.hpp"
#include "../src/integrate.hpp"
#include "../src/boundary.hpp"
#include "../src/hydrogenic.hpp"

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

double trapzIntTest(double (*f)(double), double x0 = 0, double x1 = 1, int N = 200)
{
    vector<double> x, y;

    x = linGrid(x0, x1, N);
    y = applyFunc(f, x);

    return trapzInt(x, y);
}

double shootTest(double (*fQ)(double), double (*fA)(double), double (*fB)(double),
                 double x0 = 0, double x1 = 1, int N = 200, char type = 'Q', char dir = 'f')
{
    double err = 0.0;
    vector<double> x, A(N), B(N), Q(N);

    x = linGrid(x0, x1, N);
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
        Q[N - 2] = fQ(x[N - 2]);
        Q[N - 1] = fQ(x[N - 1]);
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

    for (int i = 0; i < N; ++i)
    {
        err += abs(Q[i] - fQ(x[i]));
    }
    err /= N;

    return err;
}

double shootQPTest(double (*fQ)(double), double (*fP)(double),
                   double (*fAA)(double), double (*fAB)(double),
                   double (*fBA)(double), double (*fBB)(double),
                   double x0 = 0, double x1 = 1, int N = 200, char dir = 'f')
{
    double err = 0.0;
    vector<double> x, AA(N), AB(N), BA(N), BB(N), Q(N), P(N);

    x = linGrid(x0, x1, N);
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
        Q[N - 2] = fQ(x[N - 2]);
        Q[N - 1] = fQ(x[N - 1]);
        P[N - 2] = fP(x[N - 2]);
        P[N - 1] = fP(x[N - 1]);
    }

    shootQP(Q, P, AA, AB, BA, BB, x[1] - x[0], dir == 'f' ? -1 : 0, dir);

    for (int i = 0; i < N; ++i)
    {
        err += abs(Q[i] - fQ(x[i]));
        err += abs(P[i] - fP(x[i]));
    }
    err /= 2 * N;

    return err;
}

double shootPotentialTest(double (*fV)(double), double (*frho)(double),
                          double x0 = 0, double x1 = 1, int N = 200)
{
    double err = 0.0;
    vector<vector<double>> grid = logGrid(x0, x1, N);
    vector<double> rho(N), V(N);

    rho = applyFunc(frho, grid[1]);

    shootPotentialLog(V, rho, grid[0][1] - grid[0][0]);

    for (int i = 0; i < N; ++i)
    {
        err += abs(V[i] - fV(grid[1][i]));
    }
    err /= N;

    return err;
}

double diracTest(double Z, double mu, int n, int k,
                 double x0 = 1e-3, double x1 = 1e1, int N = 200)
{
    double err = 0.0;
    vector<double> Q(N), P(N), V(N);
    vector<vector<double>> grid = logGrid(x0, x1, N);
    double E = hydrogenicDiracEnergy(Z, mu, n, k);
    vector<vector<double>> PQ = hydrogenicDiracWavefunction(grid[1], Z, mu, n, k);

    for (int i = 0; i < N; ++i)
    {
        V[i] = -Z / grid[1][i];
    }

    boundaryDiracCoulomb(Q, P, grid[1], E, k, mu, Z);
    TurningPoint tp = shootDiracLog(Q, P, grid[1], V, E, k, mu, grid[0][1] - grid[0][0]);

    // Now on to rescale in order to make the two functions comparable
    double norm_i = tp.Pi / PQ[0][tp.i];
    double norm_e = tp.Pe / PQ[0][tp.i];

    for (int i = 0; i < N; ++i)
    {

        P[i] /= i >= tp.i ? norm_e : norm_i;
        Q[i] /= i >= tp.i ? norm_e : norm_i;

        err += abs(P[i] - PQ[0][i]);
        err += abs(Q[i] - PQ[1][i]) * Physical::c; // The c factor brings this to the same order of magnitude as P
    }
    err /= 2 * N;

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

TEST_CASE("Dirac integration", "[shootDiracLog]")
{
    // These at the moment fail, the tolerance is too low
    REQUIRE(diracTest(1, 1, 1, -1, 1e-4, 1e2, 1000) < ERRTOL);
    REQUIRE(diracTest(1, 1, 2, -1, 2e-4, 2e2, 1000) < ERRTOL);
    REQUIRE(diracTest(1, 1, 2,  1, 2e-4, 2e2, 1000) < ERRTOL);
    REQUIRE(diracTest(1, 1, 3,  1, 2e-4, 2e2, 1000) < ERRTOL);
    REQUIRE(diracTest(5, 1, 1, -1, 1e-4, 1e2, 1000) < ERRTOL);
}