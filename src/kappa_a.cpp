#include <iostream>
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // Get the atomic number
    if (argc < 1)
    {
        cout << "Must input one atomic number";
        return -1;
    }

    double Z = stod(argv[1]);
    cout << "Solving for Z = " << Z << "\n";
    DiracState s1, p2;
    DiracAtom da = DiracAtom(Z, Physical::m_mu, 2*Z*Physical::amu, 1.6);
    da.setGrid(1e-5, 1, 2000);

    try {
        s1 = da.getState(1, 0, false);
        p2 = da.getState(2, 1, false);
    }
    catch (const char *e) {
        cout << "ERROR: " << e << "\n";
    }

    cout << "Ka = " << (p2.E - s1.E) / (1e3*Physical::eV) << " keV\n";
}