#include <iostream>
#include <vector>
#include "../lib/hydrogenic.hpp"
#include "../lib/constants.hpp"
#include "../lib/utils.hpp"
#include "../vendor/easyloggingpp/easylogging++.h"

using namespace std;

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv)
{
    el::Configurations loggerConf;
    loggerConf.setToDefault();
    loggerConf.set(el::Level::Global, el::ConfigurationType::ToStandardOutput, "FALSE");
    loggerConf.set(el::Level::Global, el::ConfigurationType::Filename, "hydrogenic_print.log");
    el::Loggers::reconfigureLogger("default", loggerConf);

    LOG(INFO) << "My first info log using default logger";

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
    double r_out = 5;
    double r_in = 1e6;
    if (argc > 5)
        r_out = stod(argv[5]);
    if (argc > 6)
        r_in = stod(argv[6]);
    double m = mu ? Physical::m_mu : Physical::m_e;

    double schE = hydrogenicSchroEnergy(Z, m, n);
    double dirE = hydrogenicDiracEnergy(Z, m, n, k, true);

    double r_turn = -Z / dirE;

    double gamma = pow(k, 2) - pow(Z * Physical::alpha, 2);
    gamma = sqrt(gamma);
    double K = -dirE * Physical::alpha * (2 * m * Physical::c + dirE * Physical::alpha);
    K = sqrt(K);

    double eps = 1e-3;
    double r_low = pow(eps, 1.0 / gamma) / exp(1) * gamma / K;

    vector<vector<double>> grid = logGrid(r_low, r_turn * r_out, 1000);
    vector<double> schPsi = hydrogenicSchroWavefunction(grid[1], Z, m, n, l);
    vector<vector<double>> dirPsi = hydrogenicDiracWavefunction(grid[1], Z, m, n, k);

    for (int i = 0; i < grid[0].size(); ++i)
    {
        cout << grid[0][i] << '\t' << grid[1][i] << '\t' << schPsi[i] << '\t' << dirPsi[0][i] << '\t' << dirPsi[1][i] << '\n';
    }
}