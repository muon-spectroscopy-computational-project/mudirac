#include <iostream>
#include <vector>
#include "../src/atom.hpp"

using namespace std;

int main()
{
    Atom a = Atom();

    vector<double> r;

    r = a.getGrid();

    for (int i = 0; i < r.size(); ++i) {
        cout << r[i] << '\n';
    }
}