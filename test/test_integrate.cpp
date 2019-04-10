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
    /* 
        Q = exp(-x**2) + x
        Q' = -2x(Q-x) + 1 = (-2x)Q+(2x^2+1)
    */
    auto fQ = [](double x) { return exp(-x * x) + x; };
    auto fA = [](double x) { return -2 * x; };
    auto fB = [](double x) { return 2 * x * x + 1; };
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000) < ERRTOL);
}

TEST_CASE("Numerov integration", "[shootNumerov]")
{
    auto fQ = [](double x) { return sin(3 * x) + x * x; };
    auto fA = [](double x) { return -9.0; };
    auto fB = [](double x) { return 9.0 * x * x + 2; };
    REQUIRE(shootTest(fQ, fA, fB, 0, 1, 1000, 'N') < ERRTOL);
}

int old()
{
    int N = 1000;
    double h = 10.0 / N;
    double errQ = 0, errP = 0;

    vector<double> Q(N), P(N);
    vector<double> AA(N), AB(N), BA(N), BB(N);

    /* Basic integral test

    */

    for (int i = 0; i < N; ++i)
    {
        Q[i] = h * i + 1;
        P[i] = 1.0 / Q[i];
    }

    cout << "Error: " << (log(11.0) - trapzInt(Q, P)) / log(11.0) << '\n';

    /* Simple shootQ test 

    Q = exp(-x**2) + x
    Q' = -2x(Q-x) + 1 = -2xQ+2x^2+1

    */

    Q[0] = 1;
    Q[1] = exp(-h * h) + h;

    for (int i = 0; i < N; ++i)
    {
        AA[i] = -2 * h * i;
        AB[i] = 2 * h * h * i * i + 1;
    }

    shootQ(Q, AA, AB, h);

    // Compute error
    errQ = 0;
    for (int i = 0; i < N; ++i)
    {
        errQ += pow(Q[i] - exp(-h * h * i * i) - h * i, 2);
    }
    errQ = sqrt(errQ) / N;

    cout << "Error: Q = " << errQ << "\n";

    /*
    A simple example:

    Q' = P
    P' = -Q

    */

    Q[0] = 0;
    Q[1] = sin(h);
    P[0] = 1;
    P[1] = cos(h);

    AA = vector<double>(N, 0);
    AB = vector<double>(N, 1);
    BA = vector<double>(N, -1);
    BB = vector<double>(N, 0);

    shootQP(Q, P, AA, AB, BA, BB, h);

    // Compute error
    errQ = 0;
    errP = 0;
    for (int i = 0; i < N; ++i)
    {
        errQ += pow(Q[i] - sin(h * i), 2);
        errP += pow(P[i] - cos(h * i), 2);
        // cout << h * i << "\t" << Q[i] << "\t" << P[i] << "\n";
    }

    errQ = sqrt(errQ) / N;
    errP = sqrt(errP) / N;

    cout << "Errors: Q = " << errQ << ", P = " << errP << "\n";

    // Numerov test

    AA = vector<double>(N, 0);
    AB = vector<double>(N, 1);

    for (int i = 0; i < Q.size(); ++i)
    {
        AB[i] = i * h;
    }

    Q[0] = 0;
    Q[1] = 1.0 / 6.0 * pow(h, 3);

    shootNumerov(Q, AA, AB, h);

    // for (int i = 0; i < Q.size(); ++i)
    // {
    //     cout << h * i << '\t' << Q[i] << '\n';
    // }

    // Potential test
    vector<vector<double>> lgrid;
    vector<double> rho(N, 0), V(N, 0);

    lgrid = logGrid(1e-5, 1e1, N);

    for (int i = 0; i < N; ++i)
    {
        rho[i] *= pow(lgrid[1][i], 2);
    }
    // shootPotentialLog(V, rho, lgrid[0][1]);
    // for (int i = 0; i < N; ++i)
    // {
    //     cout << lgrid[1][i] << '\t' << V[i] << '\n';
    // }

    // vector<double> r = lgrid[1], V(r.size());

    // for (int i = 0; i < r.size(); ++i)
    // {
    //     V[i] = -1.0 / r[i];
    // }

    // double B = -0.45;
    // int k = -1;
    // double dE;

    // for (double Btest = -0.51; Btest < -0.49; Btest += 0.0001)
    // {
    //     boundaryDiracCoulomb(Q, P, r, pow(Physical::c, 2) + Btest, k);
    //     shootDiracLog(Q, P, r, V, pow(Physical::c, 2) + Btest, k, 1, lgrid[0][1] - lgrid[0][0]);
    // }

    //

    // cout << "dE = " << dE << "\n";

    // for (int i = 0; i < r.size(); ++i) {
    //     cout << r[i] << '\t' << Q[i] << '\t' << P[i] << '\n';
    // }
}