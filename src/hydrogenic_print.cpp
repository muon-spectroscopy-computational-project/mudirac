#include <iostream>
#include <vector>
#include "../lib/hydrogenic.hpp"
#include "../lib/constants.hpp"
#include "../lib/utils.hpp"

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
    int n = stoi(argv[2]);
    int k = stoi(argv[3]);
    int l = k < 0 ? -(k + 1) : k; // Angular momentum number
    bool mu = false;
    if (argc > 4)
        mu = (string(argv[4]) == "mu");
    double r_span = 5;
    if (argc > 5)
        r_span = stod(argv[5]);
    double m = mu ? Physical::m_mu : Physical::m_e;

    double schE = hydrogenicSchroEnergy(Z, m, n);
    double dirE = hydrogenicDiracEnergy(Z, m, n, k, true);

    double r_turn = -Z / dirE;

    vector<double> r = linGrid(r_turn / r_span, r_turn * r_span);
    vector<double> schPsi = hydrogenicSchroWavefunction(r, Z, m, n, l);
    vector<vector<double>> dirPsi = hydrogenicDiracWavefunction(r, Z, m, n, k);

    for (int i = 0; i < r.size(); ++i)
    {
        cout << r[i] << '\t' << schPsi[i] << '\t' << dirPsi[0][i] << '\t' << dirPsi[0][i] << '\n';
    }
}