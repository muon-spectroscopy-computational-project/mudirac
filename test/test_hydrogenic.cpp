/**
 * Test functions for hydrogenic.cpp (known solutions for hydrogen-like atoms)
 */

#include <iostream>
#include "../src/hydrogenic.hpp"
#include "../src/utils.hpp"

using namespace std;

int main() {

    double mu_H;
    double A = 16;
    int Z = 8;

    mu_H = effectiveMass(Physical::m_e, Physical::m_p*A);

    cout << "Hydrogen atom energy: " << hydrogenicSchroEnergy(Z, mu_H, 1) << "\n";
    cout << "Hydrogen atom Dirac energy: " << hydrogenicDiracEnergy(Z, mu_H, 1, -1, true)<< "\n";

    mu_H = effectiveMass(Physical::m_mu, Physical::m_p*A);

    cout << "Hydrogen muonic atom energy: " << hydrogenicSchroEnergy(Z, mu_H, 1) << "\n";
    cout << "Hydrogen muonic atom Dirac energy: " << hydrogenicDiracEnergy(Z, mu_H, 1, -1, true)<< "\n";

    // for (double r = 0.0; r < 100.0; r += 0.1) {
    //     cout << r << " " << hydrogenicSchroWavefunction(r) << " " << hydrogenicSchroWavefunction(r, 1, 1, 5, 2) << "\n";
    // }
}