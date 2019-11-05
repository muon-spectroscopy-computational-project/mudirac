#include "mudirac.hpp"

int main(int argc, char *argv[])
{
    string seed = "mudirac";
    MuDiracInputFile config;

    if (argc < 2)
    {
        cout << "Input file missing\n";
        cout << "Please use the program as `mudirac <input_file>`\n";
        cout << "Quitting...\n";
        return -1;
    }

    seed = splitString(argv[1], ".")[0];
    try
    {
        config.parseFile(argv[1]);
    }
    catch (runtime_error e)
    {
        cout << "Invalid configuration file:\n";
        cout << e.what() << "\n";
        return -1;
    }

    int output_verbosity = config.getIntValue("output");

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

    // Print out potential at high levels of verbosity
    if (output_verbosity >= 2 && (da.getPotentialFlags() && da.HAS_ELECTRONIC))
    {
        writeEConfPotential(da.getPotentialElectronic(), seed + ".epot.dat");
    }

    // Now unravel the required spectral lines
    vector<string> xr_lines = config.getStringValues("xr_lines");
    vector<int> n1range, n2range, l1range, l2range;
    vector<bool> s1range, s2range;
    vector<TransLineSpec> transqnums;

    for (int i = 0; i < xr_lines.size(); ++i)
    {
        vector<string> ranges = splitString(xr_lines[i], "-");

        if (ranges.size() != 2)
        {
            throw invalid_argument("Invalid spectral line in input file");
        }

        vector<int> nr, lr;
        vector<bool> sr;

        parseIupacRange(ranges[0], nr, lr, sr);
        n1range.insert(n1range.end(), nr.begin(), nr.end());
        l1range.insert(l1range.end(), lr.begin(), lr.end());
        s1range.insert(s1range.end(), sr.begin(), sr.end());

        parseIupacRange(ranges[1], nr, lr, sr);
        n2range.insert(n2range.end(), nr.begin(), nr.end());
        l2range.insert(l2range.end(), lr.begin(), lr.end());
        s2range.insert(s2range.end(), sr.begin(), sr.end());
    }

    for (int i = 0; i < n1range.size(); ++i)
    {
        for (int j = 0; j < n2range.size(); ++j)
        {
            TransLineSpec tnums;
            tnums.n1 = n1range[i];
            tnums.l1 = l1range[i];
            tnums.s1 = s1range[i];

            tnums.n2 = n2range[i];
            tnums.l2 = l2range[i];
            tnums.s2 = s2range[i];

            if (tnums.n2 > tnums.n1 || abs(tnums.l2 - tnums.l1) != 1)
            {
                continue;
            }

            transqnums.push_back(tnums);
        }
    }

    vector<TransitionData> transitions;

    for (int i = 0; i < transqnums.size(); ++i)
    {
        int n1, l1, n2, l2;
        bool s1, s2;
        TransitionData tdata;

        n1 = transqnums[i].n1;
        l1 = transqnums[i].l1;
        s1 = transqnums[i].s1;
        n2 = transqnums[i].n2;
        l2 = transqnums[i].l2;
        s2 = transqnums[i].s2;

        string state1, state2;
        state1 = printIupacState(n1, l1, s1);
        state2 = printIupacState(n2, l2, s2);

        LOG(INFO) << "Computing transition " << state1 << " - " << state2 << "\n";

        try
        {
            LOG(INFO) << "Computing state " << state1 << "\n";
            tdata.ds1 = da.getState(n1, l1, s1);
            LOG(INFO) << "Computing state " << state2 << "\n";
            tdata.ds2 = da.getState(n2, l2, s2);
        }
        catch (AtomErrorCode aerr)
        {
            LOG(ERROR) << SPECIAL << "Transition energy calculation for line " << state1 << "-" << state2 << " failed with AtomErrorCode " << aerr << "\n";
            return -1;
        }
        catch (const exception &e)
        {
            LOG(ERROR) << SPECIAL << "Unknown error: " << e.what() << "\n";
            return -1;
        }

        // Compute transition probability
        tdata.tmat = da.getTransitionProbabilities(n2, l2, s2, n1, l1, s1);

        LOG(INFO) << "Transition energy = " << (tdata.ds2.E - tdata.ds1.E) / (Physical::eV * 1000) << " kEv\n";

        transitions.push_back(tdata);
    }

    // Now create output files
    if (output_verbosity >= 1)
    {
        // Save a file for all lines
        ofstream out(seed + ".xr.out");

        out << "# Z = " << da.getZ() << ", A = " << da.getA() << " amu, m = " << da.getm() << " au\n";
        out << "Line\tDeltaE (eV)\tW_12 (s^-1)\n";

        for (int i = 0; i < transitions.size(); ++i)
        {
            out << xr_lines[i] << '\t' << (transitions[i].ds2.E - transitions[i].ds1.E) / Physical::eV;
            out << "\t\t" << transitions[i].tmat.totalRate() * Physical::s << '\n';
        }

        if (config.getBoolValue("write_spec"))
        {
            // Write a spectrum
            writeSimSpec(transitions, config.getDoubleValue("spec_step"), config.getDoubleValue("spec_linewidth"), config.getDoubleValue("spec_expdec"),
                         seed + ".spec.dat");
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
                string fname = seed + "." + xr_lines[i] + "." + to_string(j + 1) + ".out";
                DiracState ds = (j == 0 ? transitions[i].ds1 : transitions[i].ds2);

                LOG(DEBUG) << "Printing out state file for line " << xr_lines[i] << ", state " << (j + 1) << "\n";

                writeDiracState(ds, fname);
            }
            string fname = seed + "." + xr_lines[i] + ".tmat.out";
            writeTransitionMatrix(transitions[i].tmat, fname);
        }
    }
}