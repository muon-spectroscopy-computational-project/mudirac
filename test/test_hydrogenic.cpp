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
}