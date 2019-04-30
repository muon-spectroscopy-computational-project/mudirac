#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include "../src/atom.hpp"

using namespace std;

int main()
{
    // Atom a = Atom(1, 1, 1, 2);

    vector<double> r, V, V2;
    tuple<int, int, int> ind;
    map<tuple<int, int, int>, DiracState *> test;
    DiracAtom da = DiracAtom(1, 1, 1);
    // cout << da.getmu() << '\n';
    DiracAtom da2 = DiracAtom(1, 1, 1, 2e-2);
    DiracState s1, s1f, s2, s2f;

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
        s2 = da.getState(1, 0, false);
        s2f = da2.getState(1, 0, false);
    }
    catch (const char *s)
    {
        std::cerr << s << '\n';
    }
    // // cout << a.getZ() << '\n';
    // // cout << a.getmu() << '\n';

    for (int i = 0; i < r.size(); ++i)
    {
        cout << r[i] << '\t' << s2.P[i] << '\t' << s2f.P[i] << '\t' << s2.Q[i] << '\t' << s2f.Q[i]<< '\n';
    }

    map<int, float> testmap;

    cout << "Map test:\n";
    cout << testmap[2] << "\n"; 


}