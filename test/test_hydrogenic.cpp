/**
 * Test functions for hydrogenic.cpp (known solutions for hydrogen-like atoms)
 */

#include <iostream>
#include <vector>
#include "../src/hydrogenic.hpp"
#include "../src/utils.hpp"

using namespace std;

int main() {

    double mu;
    double A = 55.85;
    int Z = 26;
    vector<double> gf;

    mu = effectiveMass(Physical::m_e, Physical::m_p*A);

    cout << "Hydrogen atom energy: " << hydrogenicSchroEnergy(Z, mu, 1) << "\n";
    cout << "Hydrogen atom Dirac energy: " << hydrogenicDiracEnergy(Z, mu, 1, -1, true) << "\n";

    mu = effectiveMass(Physical::m_mu, Physical::m_p*A);

    cout << "Hydrogen muonic atom energy: " << hydrogenicSchroEnergy(Z, mu, 1) << "\n";
    cout << "Hydrogen muonic atom Dirac energy: " << hydrogenicDiracEnergy(Z, mu, 1, -1, true)<< "\n";

    // Try a line
    cout << "Ka: " << (hydrogenicDiracEnergy(Z, mu, 2, 1) - hydrogenicDiracEnergy(Z, mu, 1, -1))/Physical::eV*1e-3 << " keV\n";

    for (double r = 1e-5; r < 10.0; r *= 1.02) {
        gf = hydrogenicDiracWavefunction(r);
        cout << r << " " << hydrogenicSchroWavefunction(r) << " " << gf[0] << " " << gf[1] << "\n";
    }
}