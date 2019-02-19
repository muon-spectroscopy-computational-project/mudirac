#include <vector>
#include <iostream>
#include <math.h>
#include "../src/integrate.hpp"
#include "../src/boundary.hpp"

using namespace std;

int main()
{
    int N = 1000;
    double h = 10.0 / N;
    double errQ = 0, errP = 0;

    vector<double> Q(N), P(N);
    vector<double> AA(N), AB(N), BA(N), BB(N);

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
    vector<double> rho(N, 1), V(N, 0);

    lgrid = logGrid(1e-5, 1e1, N);

    for (int i = 0; i < N; ++i)
    {
        rho[i] *= pow(lgrid[1][i], 2);
    }
    shootPotentialLog(V, rho, lgrid[0][1]);
    for (int i = 0; i < N; ++i)
    {
        cout << lgrid[1][i] << '\t' << V[i] << '\n';
    }

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