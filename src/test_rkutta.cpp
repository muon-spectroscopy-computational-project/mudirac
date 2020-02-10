#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include "../lib/constants.hpp"
#include "../lib/integrate.hpp"
#include "../lib/hydrogenic.hpp"

#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

int main() {

    int N = 3000;
    int n = 1;
    int k = -1;
    double E;
    double Z = 12;
    double dx = 0.006;
    double gp, mc;
    vector<double> r(N);
    vector<vector<double>> PQ;
    vector<double> y(N), zeta_q(N), zeta_rk(N), V(N);
    vector<double> A(N, 0);
    vector<double> B(N, 0);

    r[0] = 1e-5;
    for (int i = 1; i < N; ++i) {
        r[i] = r[i-1]*exp(dx);
    }

    mc = Physical::m_mu*Physical::c;
    E = hydrogenicDiracEnergy(Z, Physical::m_mu, n, k);
    PQ = hydrogenicDiracWavefunction(r, Z, Physical::m_mu, n, k);
    
    for (int i = 0; i < N; ++i) {
        y[i] = PQ[1][i]/PQ[0][i];
        V[i] = -Z/r[i];
        gp = (mc + (E - V[i]) * Physical::alpha);
        A[i] = 2 * (k - gp * r[i] * y[i]);
        B[i] = -r[i] * (1 + y[i] * y[i]) * Physical::alpha;
    }

    shootQ(zeta_q, A, B, dx);
    shootRungeKutta(zeta_rk, A, B, dx);

    for (int i = 0; i < N; ++i) {
        cout << r[i] << "\t" << zeta_q[i] << "\t" << zeta_rk[i] << "\n";
    }

}