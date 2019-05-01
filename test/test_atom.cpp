#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <tuple>
#include <limits>
#include "../lib/atom.hpp"

using namespace std;

int main()
{
    // Atom a = Atom(1, 1, 1, 2);
    cout << log(numeric_limits<double>::min()) << "\n";
    return 0;

    vector<double> r, V, V2;
    tuple<int, int, int> ind;
    map<tuple<int, int, int>, DiracState *> test;
    DiracAtom da = DiracAtom(2, Physical::m_mu, -1);
    DiracAtom da2 = DiracAtom(2, Physical::m_mu, -1);
    DiracState s1, s1f, s2, s2f, p2, p2f;

    da.setGrid(1e-6, 1);
    da2.setGrid(1e-6, 1);

    try
    {
        da.calcState(1, 0, false);
        da.calcState(2, 0, false);
        da2.calcState(1, 0, false);
        da2.calcState(2, 0, false);
    }
    catch (const char *s)
    {
        std::cerr << s << '\n';
    }

    r = da.getGrid();
    V = da.getPotential();
    V2 = da2.getPotential();

    da.maxit = 1000;
    da2.maxit = 1000;
    try {
        s1 = da.getState(1, 0, false);
        s1f = da2.getState(1, 0, false);
        s2 = da.getState(2, 0, false);
        s2f = da2.getState(2, 0, false);
        p2 = da.getState(2, 1, false);
        p2f = da2.getState(2, 1, false);
    }
    catch (const char *s)
    {
        std::cerr << s << '\n';
    }
    // // cout << a.getZ() << '\n';
    // // cout << a.getmu() << '\n';

    cout << (s2.E - s1.E)/Physical::eV<< "\n";
    cout << (p2.E - s1.E)/Physical::eV << "\n";
    for (int i = 0; i < r.size(); ++i)
    {
        cout << r[i] << '\t' << s1.P[i] << '\t' << s2.P[i] << '\t' << p2.P[i] << '\n';
    }

}