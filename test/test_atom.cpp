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
    DiracAtom da2 = DiracAtom(1, 1, 1, 1e-2);

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

    // cout << a.getZ() << '\n';
    // cout << a.getmu() << '\n';

    // for (int i = 0; i < r.size(); ++i)
    // {
    //     cout << r[i] << '\t' << V[i] << '\t' << V2[i] << '\n';
    // }

    // DiracState *st = new DiracState();

    cout << (test[{1, 1, 1}] == NULL) << "\n";
    cout << "Finished\n";
}