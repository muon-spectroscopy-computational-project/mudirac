#include <vector>
#include <iostream>
#include <math.h>
#include "../src/integrate.hpp"

using namespace std;

int main()
{

    /*
    A simple example:

    Q' = P
    P' = -Q

    */

    int N = 1000;
    double h = 10.0/N;
    vector<double> Q(N), P(N);
    vector<double> AA(N, 0), AB(N, 1), BA(N, -1), BB(N, 0);

    Q[0] = 0;
    Q[1] = sin(h);
    P[0] = 1;
    P[1] = cos(h);

    shootQP(Q, P, AA, AB, BA, BB, h);

    // Compute error
    double errQ = 0, errP = 0;
    for (int i = 0; i < N; ++i)
    {
        errQ += pow(Q[i]-sin(h*i), 2);
        errP += pow(P[i]-cos(h*i), 2);
        // cout << h * i << "\t" << Q[i] << "\t" << P[i] << "\n";
    }

    errQ = sqrt(errQ)/N;
    errP = sqrt(errP)/N;

    // cout << "Errors: Q = " << errQ << ", P = " << errP << "\n";

    vector<vector<double>> lgrid;

    lgrid = logGrid(1e-5, 1e1, 1000);

    vector<double> r = lgrid[1], V(r.size());

    for (int i = 0; i < r.size(); ++i) {
        V[i] = -1.0/r[i];
    }

    double B = -0.5, K = sqrt(pow(B*Physical::alpha, 2)-2*B);

    Q[0] = pow(r[0], -1);
    Q[1] = pow(r[1], -1);
    P[0] = pow(r[0], 1);
    P[1] = pow(r[1], 1);

    Q[N-2] = exp(-K*r[N-2]);
    Q[N-1] = exp(-K*r[N-1]);
    P[N-2] = -Q[N-2]*(B*Physical::alpha+2*Physical::c)/K;
    P[N-1] = -Q[N-1]*(B*Physical::alpha+2*Physical::c)/K;

    int turn_i = shootDiracLog(Q, P, r, V, pow(Physical::c, 2)+B, -1, 1, lgrid[0][1]-lgrid[0][0]);

    for (int i = 0; i < r.size(); ++i) {
        if (i < turn_i) {
            Q[i] *= Q[turn_i]/Q[turn_i-1];
            P[i] *= P[turn_i]/P[turn_i-1];
        }
        cout << r[i] << '\t' << Q[i] << '\t' << P[i] << '\n';
    }

}