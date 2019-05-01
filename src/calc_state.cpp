#include <iostream>
#include <vector>
#include "../lib/atom.hpp"
#include "../lib/constants.hpp"

using namespace std;

int main(int argc, char **argv)
{
    // Parameters
    if (argc != 5) {
        cout << "Wrong number of arguments\n";
        return -1;
    }
    
    double Z = stod(argv[1]);
    int n = stoi(argv[2]);
    int l = stoi(argv[3]);
    bool s = stoi(argv[4]) != 0;
    
    DiracAtom da = DiracAtom(Z, Physical::m_mu);
    DiracState ds;

    da.setGrid(1e-5, 1);

    vector<double> r = da.getGrid();

    try {
        ds = da.getState(n, l, s);
    }
    catch (AtomConvergenceException& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return -1;
    }
    catch (const char* e) {
        cerr << "ERROR: " << e << "\n";
        return -1;
    }

    cout << "E = " << ds.E/Physical::eV << " eV\n";
    for (int i = 0; i < r.size(); ++i) {
        cout << r[i] << '\t' << ds.P[i] << '\t' << ds.Q[i] << '\n';
    }

}