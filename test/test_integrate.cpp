#include <vector>
#include <iostream>
#include <math.h>
#include "../lib/utils.hpp"
#include "../lib/integrate.hpp"
#include "../lib/boundary.hpp"
#include "../lib/hydrogenic.hpp"

#include "../vendor/catch/catch.hpp"

#define ERRTOL_VHIGH 1e-7
#define ERRTOL_HIGH 1e-5
#define ERRTOL_LOW 1e-3

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

double trapzIntTest(double (*f)(double), double x0 = 0, double x1 = 1, int N = 200, bool step = false, bool loggrid=false)
{
    double dx;
    vector<double> x, y, lx;

    double ans;
    if (loggrid) {
        vector<vector<double>> grid = logGrid(x0, x1, N);
        lx = grid[0];
        x = grid[1];
        dx = lx[1]-lx[0];
        y = applyFunc(f, x);
        y = vectorOperation(y, x, '*');
    }
    else {
        x = linGrid(x0, x1, N);
        dx = x[1]-x[0];
        y = applyFunc(f, x);
    }
    
    if (step)
    {
        ans = trapzInt(dx, y);
    }
    else
    {
        if (loggrid) {
            ans = trapzInt(lx, y);
        }
        else {
            ans = trapzInt(x, y);
        }
    }

    return ans;
}

double shootTest(double (*fQ)(double), double (*fA)(double), double (*fB)(double),
                 double x0 = 0, double x1 = 1, int N = 200, char type = 'R', char dir = 'f')
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
        if (type == 'N') // Not required for Runge-Kutta
            Q[1] = fQ(x[1]);
    }
    else
    {
        if (type == 'N')
            Q[N - 2] = fQ(x[N - 2]);
        Q[N - 1] = fQ(x[N - 1]);
    }

    switch (type)
    {
    case 'R':
        shootRungeKutta(Q, A, B, x[1] - x[0], dir == 'f' ? -1 : 0, dir);
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
        P[0] = fP(x[0]);
    }
    else
    {
        Q[N - 1] = fQ(x[N - 1]);
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

    // Apply just the known solution as boundary conditions
    for (int i = 0; i < 2; ++i)
    {
        P[i] = PQ[0][i];
        Q[i] = PQ[1][i];
        P[N - i - 1] = PQ[0][N - i - 1];
        Q[N - i - 1] = PQ[1][N - i - 1];
    }
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
    REQUIRE(trapzIntTest(cbrt, 1.0, 8.0, 200, true) == Approx(45.0 / 4.0));

    // Logarithmic integration
    REQUIRE(trapzIntTest(exp, 1e-2, 1, 1000, true, true) == Approx(exp(1) - exp(1e-2)));
}

TEST_CASE("Shooting integration", "[shootRungeKutta]")
{
    /* 
        Q = tan(x)
        Q' = tan(x)^2 + 1 = tan(x)Q+1
    */
    auto one = [](double x) { return 1.0; };
    REQUIRE(shootTest(tan, tan, one, 0, 1, 2000) < ERRTOL_VHIGH);
    REQUIRE(shootTest(tan, tan, one, 0, 1, 2000, 'R', 'b') < ERRTOL_VHIGH);
    /* 
        Q = exp(-x**2) + x
        Q' = -2x(Q-x) + 1 = (-2x)Q+(2x^2+1)
    */
    auto fQ = [](double x) { return exp(-x * x) + x; };
    auto fA = [](double x) { return -2 * x; };
    auto fB = [](double x) { return 2 * x * x + 1; };
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 2000) < ERRTOL_VHIGH);
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 2000, 'R', 'b') < ERRTOL_VHIGH);
}

TEST_CASE("Numerov integration", "[shootNumerov]")
{
    auto fQ = [](double x) { return sin(3 * x) + x * x; };
    auto fA = [](double x) { return -9.0; };
    auto fB = [](double x) { return 9.0 * x * x + 2; };
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000, 'N') < ERRTOL_HIGH);
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000, 'N', 'b') < ERRTOL_HIGH);
}

TEST_CASE("Coupled integration", "[shootQP]")
{
    auto fQ = [](double x) { return sin(x); };
    auto fP = [](double x) { return cos(x); };
    auto one = [](double x) { return 1.0; };
    auto negone = [](double x) { return -1.0; };
    auto zero = [](double x) { return 0.0; };
    REQUIRE(shootQPTest(fQ, fP, zero, one, negone, zero, 0, 1, 1000) < ERRTOL_VHIGH);
    REQUIRE(shootQPTest(fQ, fP, zero, one, negone, zero, 0, 1, 1000, 'b') < ERRTOL_VHIGH);
}

TEST_CASE("Potential integration", "[shootPotentialLog]")
{
    // Test with constant density
    auto frho = [](double x) { return 4 * M_PI * x * x; };
    auto fV = [](double x) { return 2.0 / 3.0 * M_PI * x * x; };
    REQUIRE(shootPotentialTest(fV, frho, 1e-2, 1e-2 * M_E, 1000) < ERRTOL_HIGH);
    // Test with finite sphere
    auto frho_f = [](double x) { return x < 1 ? 4 * M_PI * x * x : 0; };
    auto fV_f = [](double x) { return x < 1 ? 2.0 / 3.0 * M_PI * x * x : 2.0 / 3.0 * M_PI * (3 - 2 / x); };
    REQUIRE(shootPotentialTest(fV_f, frho_f, 1e-1, 1e1, 2000) < ERRTOL_HIGH);
}

TEST_CASE("Dirac integration", "[shootDiracLog]")
{
    // At the moment only the 'low' tolerance is achievable.
    // This needs improvement
    REQUIRE(diracTest(1, 1, 1, -1, 1e-4, 1e2, 1000) < ERRTOL_LOW);
    REQUIRE(diracTest(1, 1, 2, -1, 2e-4, 2e2, 1000) < ERRTOL_LOW);
    REQUIRE(diracTest(1, 1, 2, 1, 2e-4, 2e2, 1000) < ERRTOL_LOW);
    REQUIRE(diracTest(1, 1, 3, 1, 2e-4, 2e2, 1000) < ERRTOL_LOW);
    REQUIRE(diracTest(5, 1, 1, -1, 1e-4, 1e2, 1000) < ERRTOL_LOW);
}