#include <iostream>
#include <vector>
#include "../lib/atom.hpp"
#include "../lib/hydrogenic.hpp"
#include "../lib/constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

int main(int argc, char **argv)
{
    AixLog::Log::init({make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::normal, "calc_state.log", "#message"),
                       make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::special, "calc_state.err")});

    LOG(INFO) << "Starting calc_state...\n";

    // Parameters
    if (argc < 4)
    {
        cout << "Wrong number of arguments\n";
        return -1;
    }

    double Z = stod(argv[1]);
    int n = stoi(argv[2]);
    int l = stoi(argv[3]);
    bool s = stoi(argv[4]) != 0;
    double A = argc > 5 ? stod(argv[5]) : -1;

    int k = round((abs(l - s * 0.5) + 0.5) * (s ? 1 : -1));

    LOG(INFO) << "Calculating states for atom with Z = " << Z << ", state n = " << n << " l = " << l << " s = " << (s ? "1/2" : "-1/2") << "\n";

    vector<vector<double>> PQ;

    DiracAtom da = DiracAtom(Z, Physical::m_mu, A, NuclearRadiusModel::SPHERE);
    DiracState ds;

    // da.setGridRelative(1e-4, 1e2, 3000);

    // vector<double> r = da.getGrid();
    vector<double> r;

    try
    {
        ds = da.getState(n, l, s);
        r = ds.grid;
        cout << k << '\n';
        PQ = hydrogenicDiracWavefunction(r, Z, da.getmu(), n, k);
    }
    catch (AtomConvergenceException &e)
    {
        cerr << "ERROR: " << e.what() << "\n";
        return -1;
    }
    catch (const char *e)
    {
        cerr << "ERROR: " << e << "\n";
        return -1;
    }

    cout << "E = " << (ds.E - da.getmu() * pow(Physical::c, 2)) / Physical::eV << " eV\n";
    for (int i = 0; i < r.size(); ++i)
    {
        cout << r[i] << '\t' << ds.P[i] << '\t' << ds.Q[i] << '\t' << PQ[0][i] << '\t' << PQ[1][i] << '\n';
    }
}