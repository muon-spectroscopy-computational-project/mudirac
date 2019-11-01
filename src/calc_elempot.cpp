#include <vector>
#include <iostream>
#include "../lib/econfigs.hpp"
#include "../lib/potential.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    AixLog::Log::init<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::normal, "calc_elempot.log");

    if (argc < 2)
    {
        cout << "Type an element symbol\n";
        return -1;
    }
    string el = argv[1];
    double dx = 1e-2;
    int A = getElementMainIsotope(el);
    double R = getIsotopeRadius(el, A)*Physical::fm;

    ElectronicConfiguration econf(el, -1, 1, true, true);
    CoulombSpherePotential V_coulomb(-econf.totQ(), R);
    EConfPotential V_econf(econf, econf.innerShellRadius(), dx, 1e-4, -1, econf.outerShellRadius());
    vector<double> lgrid = V_econf.getGrid()[0];
    vector<double> grid = V_econf.getGrid()[1];
    vector<double> rho = V_econf.getrho();

    int i0 = lgrid[0]/dx;
    int i1 = lgrid.back()/dx;

    LOG(TRACE) << "Grid indices: from " << i0 << " to " << i1 << "\n";

    LOG(TRACE) << "Electronic configuration charge: " << econf.totQ() << "\n";
    LOG(TRACE) << "Total integrated charge: " << V_econf.getQ() << "\n";
    LOG(TRACE) << "Recomputed charge: " << trapzInt(dx, vectorOperation(rho, grid, '*')) << "\n";
    for (int i = i0; i <= i1; ++i)
    {
        cout << grid[i-i0] << '\t' << rho[i-i0] << '\t' << V_coulomb.V(grid[i-i0])/Physical::eV <<  '\t' << V_econf.V(grid[i-i0])/Physical::eV << '\n';
    }
}