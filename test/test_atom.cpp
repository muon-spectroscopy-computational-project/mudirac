#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include "../src/atom.hpp"

using namespace std;

int main()
{
    // Atom a = Atom(1, 1, 1, 2);

    vector<double> r, V;
    tuple<int, int, int> ind;
    map<tuple<int, int, int>, DiracState *> test;
    DiracAtom da = DiracAtom(1, 1, 1, 1e-3);

    try
    {
        da.calcState(1, 0, false);
    }
    catch (const char *s)
    {
        std::cerr << s << '\n';
    }

    // r = da.getGrid();
    // V = da.getPotential();

    // cout << a.getZ() << '\n';
    // cout << a.getmu() << '\n';

    // for (int i = 0; i < r.size(); ++i) {
    //     cout << r[i] << '\t' << V[i] << '\n';
    // }

    // DiracState *st = new DiracState();

    cout << (test[{1, 1, 1}] == NULL) << "\n";
    cout << "Finished\n";
}