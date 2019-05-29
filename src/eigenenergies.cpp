#include <iostream>
#include <fstream>
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"
#include "../lib/log.hpp"

using namespace std;

int main(int argc, char **argv)
{
    Log logger("eigenenergies.log");

    std::clog << "Test\n";

    ofstream exactE("energy_table_exact.dat"), finiteE("energy_table_finite.dat"), compareE("energy_table_1s.dat");

    double A, mu, E;
    DiracState state;
    // Test all atomic numbers, from 2 to 90, for their kappa alpha line for various models
    for (int Z = 1; Z < 70; Z += 3)
    {
        std::cout << Z << '\n';
        A = getIsotopeMass(Z);
        mu = effectiveMass(A * Physical::amu, Physical::m_mu);
        DiracAtom daCoulomb = DiracAtom(Z, Physical::m_mu, A);
        DiracAtom daSphere = DiracAtom(Z, Physical::m_mu, A, SPHERE);

        exactE << Z << '\t' << A << '\t';
        finiteE << Z << '\t' << A << '\t';
        compareE << Z << '\t' << A << '\t';

        for (int n = 1; n < 3; ++n)
        {
            try
            {
                E = hydrogenicDiracEnergy(Z, Physical::m_mu, n, -1, true);
                exactE << E / (1e3 * Physical::eV) << '\t';
                compareE << E / (1e3 * Physical::eV) << '\t';

                state = daSphere.getState(n, 0, false);
                E = state.E - mu * pow(Physical::c, 2);
                finiteE << E / (1e3 * Physical::eV) << '\t';
                compareE << E / (1e3 * Physical::eV) << '\t';

                double R = Atom::sphereNuclearModel(A);
                E = -Z/R*1.5 + 1.5*sqrt(Z/(mu*pow(R, 3)));
                compareE << E / (1e3 * Physical::eV) << '\t';
            }
            catch (TurningPointError tpe)
            {
                cout << tpe.what() << '\n';
                throw("Failed");
            }
        }

        exactE << '\n';
        finiteE << '\n';
        compareE << '\n';
    }
}