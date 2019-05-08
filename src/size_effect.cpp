#include <iostream>
#include <vector>
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // Parameters
    cout << argc << '\n';
    if (argc < 3) {
        cout << "Wrong number of arguments\n";
        return -1;
    }
    
    double Z = stod(argv[1]);
    int maxn = stoi(argv[2]);
    double A = getIsotopeMass(Z);
    double restE = Physical::m_mu*pow(Physical::c, 2.0);
    DiracState s, s_sphere;

    DiracAtom da = DiracAtom(Z, Physical::m_mu, A, NuclearRadiusModel::POINT);
    DiracAtom da_sphere = DiracAtom(Z, Physical::m_mu, A, NuclearRadiusModel::SPHERE);
    
    da.setGridRelative(1e-2, 1e2, 3000);
    da_sphere.setGridRelative(1e-2, 1e2, 3000);

    for (int n = 1; n <= maxn; ++n) {
        try {
            s = da.getState(n, 0, false);
            s_sphere = da_sphere.getState(n, 0, false);
        }
        catch (TurningPointError tpe) {
            cerr << tpe.what() << '\n';
            continue;
        }

        cout << n << '\t' << (s.E-restE)/Physical::eV << '\t' << (s_sphere.E-restE)/Physical::eV << '\n';
    }

}