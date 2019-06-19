#include <iostream>
#include <vector>
#include "../lib/utils.hpp"
#include "../lib/potential.hpp"
#include "../lib/hydrogenic.hpp"
#include "../lib/integrate.hpp"
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"

using namespace std;

int main(int argc, char **argv)
{

    // Parameters
    if (argc < 3)
    {
        std::cout << "Wrong number of arguments\n";
        return -1;
    }

    double Z = stod(argv[1]);
    int maxn = stoi(argv[2]);
    double A = getIsotopeMass(Z);
    double restE = Physical::m_mu * pow(Physical::c, 2.0);
    double R, corrE;
    vector<vector<double>> PQ, n_grid;
    CoulombSpherePotential Vpoint, Vsphere;
    vector<double> rho, dV;
    DiracState s, s_sphere;

    DiracAtom da = DiracAtom(Z, Physical::m_mu, A, NuclearRadiusModel::POINT);
    DiracAtom da_sphere = DiracAtom(Z, Physical::m_mu, A, NuclearRadiusModel::SPHERE);

    R = da_sphere.getR();
    Vpoint = CoulombSpherePotential(Z);
    Vsphere = CoulombSpherePotential(Z, R);

    for (int n = 1; n <= maxn; ++n)
    {

        try
        {
            s = da.getState(n, 0, false);
            s_sphere = da_sphere.getState(n, 0, false);
        }
        catch (TurningPointError tpe)
        {
            cerr << tpe.what() << '\n';
            continue;
        }

        n_grid = logGrid(R * 1e-4, R, 1000);

        PQ = hydrogenicDiracWavefunction(n_grid[0], Z, Physical::m_mu, n, -1);
        rho = vectorOperation(vectorOperation(PQ[0], 2.0, '^'), vectorOperation(PQ[1], 2.0, '^'), '+');

        dV = vector<double>(rho.size());
        for (int i = 0; i < dV.size(); ++i)
        {
            dV[i] = (Vsphere.V(n_grid[1][i]) - Vpoint.V(n_grid[1][i])) * n_grid[1][i];
        }

        // cout << trapzInt(n_grid[0], vectorOperation(rho, n_grid[1], '*')) << '\n';
        corrE = trapzInt(n_grid[0], vectorOperation(rho, dV, '*'));

        // PQ = hydrogenicDiracWavefunction(s_sphere.grid, Z, Physical::m_mu, n, -1);
        // rho = vectorOperation(vectorOperation(PQ[0], 2.0, '^'), vectorOperation(PQ[1], 2.0, '^'), '+');

        // cout << trapzInt(s_sphere.loggrid, vectorOperation(rho, s_sphere.grid, '*')) << '\n';
        std::cout << n << '\t' << (s.E - restE) / Physical::eV << '\t' << (s_sphere.E - restE) / Physical::eV << '\t' << corrE/Physical::eV << '\n';
    }
}