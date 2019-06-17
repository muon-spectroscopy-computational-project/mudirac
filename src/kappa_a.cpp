#include <iostream>
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

int main(int argc, char **argv)
{
    double A, mu, Ka_schro, Ka_dirac_hydro, Ka_dirac_hydro_mu, Ka_dirac_num_mu, Ka_dirac_num_sphere;
    DiracState p2, s1;

    AixLog::Log::init({make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::normal, "kappa_a.log", "#message"),
                       make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::special, "kappa_a.err")});

    // Test all atomic numbers, from 2 to 90, for their kappa alpha line for various models
    for (int Z = 1; Z < 90; ++Z)
    {
        A = getIsotopeMass(Z);
        mu = effectiveMass(A * Physical::amu, Physical::m_mu);
        DiracAtom daCoulomb = DiracAtom(Z, Physical::m_mu, A);
        DiracAtom daSphere = DiracAtom(Z, Physical::m_mu, A, SPHERE);

        try
        {
            Ka_schro = hydrogenicSchroEnergy(Z, Physical::m_mu, 2) - hydrogenicSchroEnergy(Z, Physical::m_mu, 1);
            Ka_dirac_hydro = hydrogenicDiracEnergy(Z, Physical::m_mu, 2, 1, true) - hydrogenicDiracEnergy(Z, Physical::m_mu, 1, -1, true);
            Ka_dirac_hydro_mu = hydrogenicDiracEnergy(Z, mu, 2, 1, true) - hydrogenicDiracEnergy(Z, mu, 1, -1, true);

            std::clog << "Exact energies computed\n";
            std::clog.flush();

            p2 = daCoulomb.getState(2, 1, false);
            s1 = daCoulomb.getState(1, 0, false);

            Ka_dirac_num_mu = (p2.E - s1.E);

            std::clog << "Point energies computed\n";
            std::clog.flush();

            p2 = daSphere.getState(2, 1, false);
            s1 = daSphere.getState(1, 0, false);

            Ka_dirac_num_sphere = (p2.E - s1.E);

            std::clog << "Finite nucleus energies computed\n";
            std::clog.flush();
        }
        catch (TurningPointError tpe)
        {
            cout << tpe.what() << '\n';
            throw("Failed");
        }

        // Convert to keV
        Ka_schro /= 1e3 * Physical::eV;
        Ka_dirac_hydro /= 1e3 * Physical::eV;
        Ka_dirac_hydro_mu /= 1e3 * Physical::eV;
        Ka_dirac_num_mu /= 1e3 * Physical::eV;
        Ka_dirac_num_sphere /= 1e3 * Physical::eV;

        cout << Z << '\t' << A << '\t' << Ka_schro << '\t' << Ka_dirac_hydro << '\t' << Ka_dirac_hydro_mu << '\t' << Ka_dirac_num_mu << '\t' << Ka_dirac_num_sphere << '\n';
    }
}