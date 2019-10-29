#include <vector>
#include <iostream>
#include "../lib/econfigs.hpp"
#include "../lib/potential.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    AixLog::Log::init<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::normal, "print_econf.log");

    if (argc < 2)
    {
        cout << "Insert an electronic configuration string\n";
        return -1;
    }
    string ecstr = argv[1];
    double dx = 1e-3;
    if (argc >= 3)
        dx = stod(argv[2]);

    ElectronicConfiguration econf(ecstr, -1, 1, true);
    EConfPotential V_econf(econf, 1.0, dx, 1e-3);
    vector<double> grid = V_econf.getGrid()[1];
    vector<double> rho = V_econf.getrho();

    LOG(TRACE) << "Total integrated charge: " << V_econf.getQ() << "\n";
    for (int i = 0; i < grid.size(); ++i)
    {
        cout << grid[i] << '\t' << rho[i] << '\t' << V_econf.V(grid[i]) << '\n';
    }
}