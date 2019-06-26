#include "mudirac.hpp"

int main(int argc, char *argv[])
{
    string seed = "mudirac";
    MuDiracInputFile config;

    if (argc > 1)
    {
        config.parseFile(argv[1]);
        seed = splitString(argv[1], ".")[0];
    }

    // Set up logging
    AixLog::Severity log_verbosity;
    switch (config.getIntValue("verbosity"))
    {
    case 1:
        log_verbosity = AixLog::Severity::info;
        break;
    case 2:
        log_verbosity = AixLog::Severity::debug;
        break;
    case 3:
        log_verbosity = AixLog::Severity::trace;
        break;
    default:
        log_verbosity = AixLog::Severity::info;
        break;
    }
    AixLog::Log::init({make_shared<AixLog::SinkFile>(log_verbosity, AixLog::Type::normal, seed + ".log"),
                       make_shared<AixLog::SinkFile>(AixLog::Severity::warning, AixLog::Type::special, seed + ".err")});

    LOG(INFO) << "MuDirac, a muonic atomic solver\n";
    LOG(INFO) << "by Simone Sturniolo\n";
    LOG(INFO) << "Released under the MIT License (2019)\n";
    LOG(INFO) << " \n";

    DiracAtom da = config.makeAtom();

    // Now unravel the required spectral lines
    vector<string> xr_lines = config.getStringValues("xr_lines");
    vector<pair<DiracState, DiracState>> trans_states;

    for (int i = 0; i < xr_lines.size(); ++i)
    {
        int n1, l1, n2, l2;
        bool s1, s2;
        DiracState ds1, ds2;
        vector<string> states = splitString(xr_lines[i], "-");

        if (states.size() != 2)
        {
            throw invalid_argument("Invalid spectral line in input file");
        }

        parseIupacState(states[0], n1, l1, s1);
        parseIupacState(states[1], n2, l2, s2);

        LOG(INFO) << "Computing transition " << states[0] << " - " << states[1] << "\n";

        ds1 = da.getState(n1, l1, s1);
        ds2 = da.getState(n2, l2, s2);

        LOG(INFO) << "Transition energy = " << (ds2.E - ds1.E) / (Physical::eV * 1000) << " kEv\n";

        trans_states.push_back(make_pair(ds1, ds2));
    }

    int output_verbosity = config.getIntValue("output");

    // Now create output files
    if (output_verbosity >= 1)
    {
        // Save a file for all lines
        ofstream out(seed + ".xr.out");
        DiracState ds1, ds2;

        out << "# Z = " << da.getZ() << ", A = " << da.getA() << " amu, m = " << da.getm() << " au\n";

        for (int i = 0; i < xr_lines.size(); ++i)
        {
            ds1 = trans_states[i].first;
            ds2 = trans_states[i].second;
            out << xr_lines[i] << '\t' << (ds2.E - ds1.E) / Physical::eV << '\n';
        }

        out.close();
    }

    if (output_verbosity >= 2)
    {
        // Save each individual state
        for (int i = 0; i < xr_lines.size(); ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                ofstream out(seed + "." + xr_lines[i] + "." + to_string(j + 1) + ".out");
                DiracState ds = (j == 0 ? trans_states[i].first : trans_states[i].second);

                LOG(DEBUG) << "Printing out state file for line " << xr_lines[i] << ", state " << (j+1) << "\n";

                out << "# " << ds.E / Physical::eV << '\n';
                out << "# " << ds.getn() << '\t' << ds.getl() << '\t' << ds.gets() << '\n';
                out << "# " << ds.nodes << '\t' << ds.nodesQ << '\n';

                for (int k = 0; k < ds.grid.size(); ++k)
                {
                    out << ds.grid[k] << '\t' << ds.V[k] / Physical::eV << '\t' << ds.P[k] << '\t' << ds.Q[k] << '\n';
                }

                out.close();
            }
        }
    }
}