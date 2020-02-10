#include "mudirac.hpp"

int main(int argc, char *argv[])
{
    MuDiracInputFile config;
    double E, dE;
    int k;
    DiracState ds;
    DiracAtom da;
    TurningPoint tp;
    string state_name;

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

    // Ask for further data
    cout << "Binding energy?\n";
    cin >> E;
    cout << "Quantum number k?\n";
    cin >> k;


    da = config.makeAtom();
    da.write_debug = true;

    ds = da.initState(E+da.getRestE(), k);
    da.integrateState(ds, tp, dE);
    ds.continuify(tp);
    ds.normalize();

    cout << "Turning point found at " << ds.grid[tp.i] << " atomic units\n";
    cout << "Energy error " << dE << "\n";
    state_name =  printIupacState(ds.getn(), ds.getl(), ds.gets());

    writeDiracState(ds, state_name + "_PQ.dat");
}