#include <iostream>
#include <vector>
#include "../lib/atom.hpp"
#include "../lib/hydrogenic.hpp"
#include "../lib/constants.hpp"

using namespace std;

int main(int argc, char **argv)
{

    // Parameters
    if (argc < 4)
    {
        cout << "Wrong number of arguments\n";
        return -1;
    }

    double Z = stod(argv[1]);
    int l = stoi(argv[2]);
    bool s = stoi(argv[3]) != 0;
    double A = argc > 4 ? stod(argv[4]) : -1;

    int k = round((abs(l - s * 0.5) + 0.5) * (s ? 1 : -1));

    DiracAtom da = DiracAtom(Z, Physical::m_mu, A, NuclearRadiusModel::SPHERE);
    DiracState ds;
    TurningPoint tp;
    pair<int, int> glim;

    ds.k = k;

    double E0 = hydrogenicDiracEnergy(Z, da.getmu(), 1, -1, true);
    // std::cout << E0 << '\t' << da.getV(0) << '\n';

    // da.searchBasinE(-1, 2, E0 * 2, 0);
    return 0;

    for (double fB = 5.0; fB > 1e-4; fB -= 5e-3)
    {
        try
        {
            double B = E0 * fB;
            if (B < da.getV(0))
                continue;
            double E = B + da.getRestE();
            glim = da.gridLimits(E, k);
            ds = DiracState(da.getrc(), da.getdx(), glim.first, glim.second);
            ds.k = k;
            ds.E = E;
            ds.V = da.getV(ds.grid);
            double dE;
            // dE = da.stateIntegrate(ds, tp);

            std::cout << B << '\t' << dE << '\n';
        }
        catch (TurningPointError tpe)
        {
            std::cerr << tpe.what() << '\n';
            return 1;
        }
        catch (char const *e)
        {
            std::cerr << e << '\n';
            return 1;
        }
    }
}