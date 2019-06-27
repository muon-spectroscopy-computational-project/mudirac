#include <iostream>
#include <vector>
#include "../lib/atom.hpp"
#include "../lib/config.hpp"
#include "../lib/output.hpp"
#include "../lib/hydrogenic.hpp"
#include "../lib/constants.hpp"
#include "../vendor/aixlog/aixlog.hpp"

using namespace std;

int main(int argc, char **argv)
{
    AixLog::Log::init({make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::normal, "arbitrary_state.log"),
                       make_shared<AixLog::SinkFile>(AixLog::Severity::trace, AixLog::Type::special, "arbitrary_state.err")});

    LOG(INFO) << "Starting arbitrary_state...\n";

    // Parameters
    MuDiracInputFile config;
    config.parseFile(argv[1]);
    int l = stoi(argv[2]);
    bool s = stoi(argv[3]) != 0;
    double E = stod(argv[4]);

    int k;

    DiracAtom da = config.makeAtom();
    E += da.getRestE();
    qnumSchro2Dirac(l, s, k);

    cout << hydrogenicDiracEnergy(da.getZ(), da.getmu(), 2, k, true) << "\n";

    DiracState ds = da.initState(E, k);
    TurningPoint tp;
    try
    {
        da.integrateState(ds, tp);
    }
    catch (const char *err)
    {
        LOG(ERROR) << err << "\n";
    }
    ds.continuify(tp);
    ds.normalize();
    ds.findNodes(config.getDoubleValue("node_tol"));
    LOG(TRACE) << "State converged, nodes = " << ds.nodes << ", nodesQ = " << ds.nodesQ << "\n";

    writeDiracState(ds, "arbitrary.dat");
}