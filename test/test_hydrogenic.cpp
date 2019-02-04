/**
 * Test functions for hydrogenic.cpp (known solutions for hydrogen-like atoms)
 */

#include <iostream>
#include "../src/hydrogenic.hpp"
#include "../src/utils.hpp"

using namespace std;

int main() {

    double mu_H;

    mu_H = effectiveMass(Physical::m_e, Physical::m_p);

    cout << "Hydrogen atom energy: " << hydrogenicSchroEnergy(1, mu_H, 1) << "\n";

    mu_H = effectiveMass(Physical::m_mu, Physical::m_p);

    cout << "Hydrogen muonic atom energy: " << hydrogenicSchroEnergy(1, mu_H, 1) << "\n";

    for (double r = 0.0; r < 100.0; r += 0.1) {
        cout << r << " " << hydrogenicSchroWavefunction(r) << " " << hydrogenicSchroWavefunction(r, 1, 1, 5, 2) << "\n";
    }
}