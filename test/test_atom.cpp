#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include "../src/atom.hpp"

using namespace std;

int main()
{
    Atom a = Atom(1, 1, 1, 2);

    vector<double> r, V;
    tuple<int, int, int> ind;
    map<tuple<int, int, int>, vector<double>> test;

    test[{1,1,1}] = V;

    cout << (test[{3,3,1}].empty()) << '\n'; 

    r = a.getGrid();
    V = a.getPotential();

    // cout << a.getZ() << '\n';
    // cout << a.getmu() << '\n';

    // for (int i = 0; i < r.size(); ++i) {
    //     cout << r[i] << '\t' << V[i] << '\n';
    // }
}