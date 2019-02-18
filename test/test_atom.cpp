#include <iostream>
#include <vector>
#include "../src/atom.hpp"

using namespace std;

int main()
{
    Atom a = Atom(1, 1, 1, 2);

    vector<double> r, V;

    r = a.getGrid();
    V = a.getPotential();

    // cout << a.getZ() << '\n';
    // cout << a.getmu() << '\n';

    for (int i = 0; i < r.size(); ++i) {
        cout << r[i] << '\t' << V[i] << '\n';
    }
}